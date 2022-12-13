#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Based on:
https://kubernetes.io/docs/tasks/job/fine-parallel-processing-work-queue/
http://peter-hoffmann.com/2012/python-simple-queue-redis-queue.html
http://redis.io/commands/rpoplpush
https://git.lsd.ufcg.edu.br/issre-tutorial/radiomics/-/blob/radiomics-app/app/rediswq.py
"""

import json
import uuid
import redis
import hashlib
from time import sleep
from multiprocessing import cpu_count


class RedisMessageQueue:
    """Simple Finite Work Queue with Redis Backend

    This work queue is finite: as long as no more work is added
    after workers start, the workers can detect when the queue
    is completely empty.

    The items in the work queue are assumed to have unique values.

    This object is not intended to be used by multiple threads
    concurrently.
    """
    def __init__(self, name: str, **redis_kwargs):
        """The default connection parameters are: host='localhost', port=6379, db=0

        The work queue is identified by "name".  The library may create other
        keys with "name" as a prefix.
        """
        self._db = redis.StrictRedis(**redis_kwargs)
        # The session ID will uniquely identify this "worker".
        self._session = str(uuid.uuid4())
        # Work queue is implemented as two queues: main, and processing.
        # Work is initially in main, and moved to processing when a client picks it up.
        self._main_q_key = name
        self._processing_q_key = name + ":processing"
        self._lease_key_prefix = name + ":leased_by_session:"
        self._result_q_key = name + ":results"
        print("Connected to Redis DB queue: {}".format(self._main_q_key))

    def session_id(self):
        """Return the ID for this session."""
        return self._session

    def _main_qsize(self):
        """Return the size of the main queue."""
        return self._db.llen(self._main_q_key)

    def _processing_qsize(self):
        """Return the size of the main queue."""
        return self._db.llen(self._processing_q_key)

    def empty(self):
        """Return True if the queue is empty, including work being done, False otherwise.

        False does not necessarily mean that there is work available to work on right now,
        """
        return self._main_qsize() == 0 and self._processing_qsize() == 0

    def check_expired_leases(self):
        """Return items whose lease expired to the main queue."""
        processing = self._db.lrange(self._processing_q_key, 0, -1)
        for item in processing:
            if not self._lease_exists(item):
                # atomically return the item to the main queue.
                # TODO(clenimar): watch leases and queue changes
                with self._db.pipeline() as pipe:
                    pipe.rpush(self._main_q_key, item)
                    pipe.lrem(self._processing_q_key, 0, item)
                    pipe.execute()

    def publish_result(self, value):
        return self._db.rpush(self._result_q_key, value)

    @staticmethod
    def _itemkey(item):
        """Returns a string that uniquely identifies an item (bytes)."""
        return hashlib.sha224(item).hexdigest()

    def _lease_exists(self, item):
        """True if a lease on 'item' exists."""
        return self._db.exists(self._lease_key_prefix + self._itemkey(item))

    def lease(self, lease_secs=60, is_blocking=True, timeout=None):
        """Begin working on an item the work queue.

        Lease the item for lease_secs.  After that time, other
        workers may consider this client to have crashed or stalled
        and pick up the item instead.

        If optional args block is true and timeout is None (the default), block
        if necessary until an item is available."""
        if is_blocking:
            item = self._db.blmove(
                first_list=self._main_q_key,
                second_list=self._processing_q_key,
                timeout=timeout,
                src="LEFT",
                dest="RIGHT",
            )
        else:
            item = self._db.lmove(
                first_list=self._main_q_key,
                second_list=self._processing_q_key,
                src="LEFT",
                dest="RIGHT",
            )
        if item:
            # Record that we (this session id) are working on a key.  Expire that
            # note after the lease timeout.
            # Note: if we crash at this line of the program, then GC will see no lease
            # for this item a later return it to the main queue.
            itemkey = self._itemkey(item)
            self._db.setex(self._lease_key_prefix + itemkey, lease_secs, self._session)
        return item

    def complete(self, value):
        """Complete working on the item with 'value'.

        If the lease expired, the item may not have completed, and some
        other worker may have picked it up.  There is no indication
        of what happened.
        """
        self._db.lrem(self._processing_q_key, 0, value)
        # If we crash here, then the GC code will try to move the value, but it will
        # not be here, which is fine.  So this does not need to be a transaction.
        itemkey = self._itemkey(value)
        self._db.delete(self._lease_key_prefix + itemkey)

    def flush(self, target: str = "db"):
        out = ""
        if target == "db":
            out = self._db.flushdb()
        elif target == "all":
            out = self._db.flushall()
        else:
            print("Unknown target: '{}'".format(target))
        print("Flushed Redis DB: {}".format(out))

    def unleash(self, value):
        """
        Unlease item becuse errors occured.
        """
        self._db.lrem(self._processing_q_key, 0, value)
        self._db.rpush(self._main_q_key, value)
        itemkey = self._itemkey(value)
        self._db.delete(self._lease_key_prefix + itemkey, self._session)
        pass

    def push(self, value: str):
        """Sloppily insert an item in the main queue."""
        o = self._db.rpush(self._main_q_key, value)
        print(f"Item '{value}' pushed to the queue '{self._main_q_key}' with number {o}")

    def insert_queue(self, value: str, queue: str):
        """Sloppily insert an item in the main queue."""
        self._db.rpush(value, queue)

    def disconnect(self):
        self._db.connection_pool.disconnect()


def strings_to_redis(
        strings: list,
        queue_name: str = "redis",
        is_flush: bool = False,
        **redis_kwargs
):
    mq = RedisMessageQueue(name=queue_name, **redis_kwargs)
    if is_flush:
        mq.flush("db")
    c = 0
    for s in strings:
        mq.push(s)
        c += 1
    print("Uploaded {} items to the Redis queue '{}' ".format(c, queue_name))
    mq.disconnect()


def redis_to_strings(
        queue_name: str = "redis",
        content_length: int = cpu_count(),
        pause: int = 10,
        **redis_kwargs
):
    mq = RedisMessageQueue(name=queue_name, **redis_kwargs)
    max_idle_counter = 10
    c = 0
    out = list()
    while c < max_idle_counter and len(out) < content_length:
        if mq.empty():
            print("The queue '{}' is empty, paused for {} seconds, {} attempts left".format(
                queue_name, pause, max_idle_counter - c)
            )
            mq.disconnect()
            mq = RedisMessageQueue(name=queue_name, **redis_kwargs)
            c += 1
            sleep(60)
        else:
            c = 0
            q = mq.lease(lease_secs=10, is_blocking=True, timeout=2)
            mq.complete(q)
            out.append(q.decode("utf-8"))
        sleep(pause)
    mq.disconnect()
    print("Fetched {} items from the Redis queue '{}' ".format(len(out), queue_name))
    return out


def dicts_to_redis(dicts: list, **kwargs):
    from meta.utils.primitive import dicts_to_strings
    strings_to_redis(strings=dicts_to_strings(dicts), **kwargs)


def redis_to_dicts(**kwargs):
    from meta.utils.primitive import strings_to_dicts
    out = redis_to_strings(**kwargs)
    return strings_to_dicts(out)
