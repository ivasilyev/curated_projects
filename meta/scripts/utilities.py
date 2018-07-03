#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def ends_with_slash(string):
    return (string + "/", string)[string.endswith("/")]


def get_page(url):
    import requests
    header = "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US) AppleWebKit/525.19 (KHTML, like Gecko) Chrome/1.0.154.53 Safari/525.19"
    return requests.get(url, headers={'User-Agent': header}).content


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


def dict2pd_series(dictionary):
    import pandas as pd
    output = pd.Series()
    for key in dictionary:
        output.at[key] = dictionary[key]
    return output


def multi_core_queue(func, queue):
    import multiprocessing
    pool = multiprocessing.Pool()
    output = pool.map(func, queue)
    pool.close()
    pool.join()
    return output


def single_core_queue(func, queue):
    return [func(i) for i in queue]


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


def get_script_dir():
    import os
    import sys
    return os.path.dirname(os.path.realpath(sys.argv[0]))
