#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os


def load_string(file: str):
    with open(file=file, mode="r", encoding="utf-8") as f:
        s = f.read()
        f.close()
    return s


def dump_string(string: str, file: str, append: bool = False):
    os.makedirs(os.path.dirname(file), exist_ok=True)
    mode = "w"
    if append:
        mode = "a"
    with open(file=file, mode=mode, encoding="utf-8") as f:
        f.write(string)
        f.close()


def load_list(file: str):
    from meta.utils.primitive import split_lines
    return split_lines(load_string(file))


def dump_list(lst: list, file: str):
    dump_string(string="\n".join([str(i) for i in lst]) + "\n", file=file)


def load_2d_array(file: str):
    from meta.utils.primitive import string_to_2d_array
    return string_to_2d_array(load_string(file))


def dump_2d_array(array: list, file: str):
    dump_list(lst=["\t".join([str(j) for j in i]) for i in array], file=file)


def load_dict(file: str):
    import json
    return json.loads(load_string(file))


def dump_dict(d: dict, file: str, **kwargs):
    _kwargs = dict(indent=4, sort_keys=False)
    if len(kwargs.keys()) > 0:
        _kwargs.update(kwargs)
    import json
    return dump_string(json.dumps(d, **_kwargs), file)
