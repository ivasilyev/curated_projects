
from datetime import datetime


def get_timestamp(x=None, fmt: str = "%Y-%m-%d-%H-%M-%S"):
    if x is None:
        return datetime.now().strftime(fmt)
    return datetime.strftime(x, fmt)


def count_elapsed_seconds(t):
    from time import perf_counter
    return f"{perf_counter() - t :.3f} s."
