__name__ = "meta"


def get_logger():
    from logging import basicConfig, getLogger, DEBUG
    basicConfig(format="%(asctime)s %(message)s")
    logger = getLogger()
    logger.setLevel(DEBUG)
    return logger
