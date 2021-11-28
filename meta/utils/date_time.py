
def get_timestamp(x=None, fmt: str = "%Y-%m-%d-%H-%M-%S"):
    from datetime import datetime
    if x is None:
        return datetime.now().strftime(fmt)
    return datetime.strftime(x, fmt)
