#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re


def remove_empty_values(input_list):
    output_list = []
    if input_list is not None:
        for i in input_list:
            if i is not None:
                try:
                    if len(i) > 0:
                        output_list.append(i)
                except TypeError:
                    continue
    return output_list


def safe_findall(pattern, string, idx: int = 0, verbose: bool = False):
    out = ""
    try:
        out = re.findall(pattern, string)[idx]
        if not isinstance(out, str):
            print(pattern, string, out)
    except IndexError:
        if verbose:
            print(f"Can't find the regex pattern '{pattern}' within the string: '{string}'")
    return out


def flatten_2d_array(array: list):
    return [j for i in array for j in i]


def split_lines(string: str):
    out = [i.strip() for i in re.split("[\r\n]+", string)]
    return remove_empty_values(out)


def string_to_2d_array(string: str):
    out = [[j.strip() for j in i.split("\t")] for i in split_lines(string)]
    return remove_empty_values(out)


def object_to_dict(o):
    out = dict()
    for key, value in o.__dict__.items():
        class_name = value.__class__.__name__
        if key.startswith("_") or class_name in ["function", "getset_descriptor"]:
            continue
        if class_name in ["property", ]:
            out[key] = o.__dict__[key].__get__(o)
            continue
        out[key] = value
    return out
