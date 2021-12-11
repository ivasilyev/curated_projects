#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def single_core_queue(func, queue) -> list:
    return [func(i) for i in queue]


def multi_core_queue(func, queue: list, processes: int = 0, async_: bool = False) -> list:
    import multiprocessing as mp
    if processes == 0:
        processes = mp.cpu_count()
    pool = mp.Pool(processes=processes)
    if async_:
        result = pool.map_async(func, queue)
    else:
        result = pool.map(func, queue)
    pool.close()
    pool.join()
    if async_:
        return result.get()
    return result


def wrapper(d: dict):
    """
    An ultimate wrapper
    :param d: Dictionary {'func': a function or a method, 'args': list, 'kwargs': dict}
    :return: The result of evaluating func(*args, **kwargs)
    """
    if "args" not in d.keys():
        d["args"] = []
    if "kwargs" not in d.keys():
        d["kwargs"] = dict()
    return d["func"](*d["args"], **d["kwargs"])


def attempt_func(func, args, raised: bool = False):
    from meta.utils.date_time import randomize_sleep
    _ATTEMPTS = 5
    attempt = 1
    while attempt <= _ATTEMPTS:
        try:
            if any(isinstance(args, i) for i in (list, tuple)):
                return func(*args)
            if any(isinstance(args, i) for i in (dict,)):
                return func(**args)
        except Exception as e:
            print("Caught exception for attempt {}: `{}`".format(attempt, e))
            attempt += 1
            if raised:
                raise
            randomize_sleep()
    print("Exceeded number of attempts for the function: '{}'".format(func.__name__))
