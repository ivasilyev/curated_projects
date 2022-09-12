
from datetime import datetime


def get_timestamp(x=None, fmt: str = "%Y-%m-%d-%H-%M-%S"):
    if x is None:
        return datetime.now().strftime(fmt)
    return datetime.strftime(x, fmt)


def count_elapsed_seconds(t) -> str:
    """
    :param t: Output of perf_counter()
    :return: Formatted string
    """
    from time import perf_counter
    return f"{perf_counter() - t :.3f} s."


def randomize_sleep(min_: int = 30, max_: int = 120):
    from time import sleep
    from random import randint
    sleep(randint(min_, max_))
