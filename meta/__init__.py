__name__ = "meta"

LOGGING_FORMAT = "%(asctime)s.%(msecs)06d [%(filename)s %(name)s] %(levelname)s %(funcName)s(%(lineno)s) - %(message)s"
DATE_FORMAT = "%Y.%m.%d %H:%M:%S"
LOGGING_LEVELS = tuple([i * 10 for i in range(6)])


def get_logger(name: str, level: int):
    # _LOG = get_logger()
    import logging as g
    from sys import stdout
    assert level in LOGGING_LEVELS
    logger = g.getLogger(name)
    logger.setLevel(level)
    formatter = g.Formatter()
    stdout_handler = g.StreamHandler(stdout)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)
    return logger
