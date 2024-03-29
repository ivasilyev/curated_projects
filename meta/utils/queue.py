#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from joblib import delayed, Parallel
from multiprocessing import cpu_count, Pool
from typing import Callable, Iterable


def single_core_queue(func: Callable, queue: Iterable) -> list:
    return [func(i) for i in queue]


def multi_core_queue(func: Callable, queue: Iterable, processes: int = 0, async_: bool = False) -> list:
    if processes == 0:
        processes = cpu_count()
    pool = Pool(processes=processes)
    if async_:
        result = pool.map_async(func, queue)
    else:
        result = pool.map(func, queue)
    pool.close()
    pool.join()
    if async_:
        return result.get()
    return result


def multi_core_queue2(func: Callable, queue: Iterable, **kwargs):
    out = Parallel(n_jobs=-1)(delayed(func)(i, **kwargs) for i in queue)
    return out


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


def attempt_func(func, attempts: int = 5, raised: bool = False, *args, **kwargs):
    for _attempt in range(attempts):
        attempt = _attempt + 1
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print("Caught exception for attempt {} of {}: `{}`".format(attempt, attempts, e))
            if raised:
                raise
    print("Exceeded number of attempts for the function: '{}'".format(func.__name__))
