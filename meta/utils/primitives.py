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


def safe_findall(pattern, string, idx: int = 0):
    try:
        return re.findall(pattern, string)[idx]
    except IndexError:
        print("Can't find the regex pattern '{}' within the string: '{}'".format(pattern, string))
        return ""


def flatten_2d_array(array: list):
    return [j for i in array for j in i]


def split_lines(string: str):
    from re import split
    out = [i.strip() for i in split("[\r\n]+", string)]
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
