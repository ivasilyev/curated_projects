#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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


def split_lines(string: str):
    from re import split
    out = [i.strip() for i in split("[\r\n]+", string)]
    return remove_empty_values(out)


def string_to_2d_array(string: str):
    out = [[j.strip() for j in i.split("\t")] for i in split_lines(string)]
    return remove_empty_values(out)

