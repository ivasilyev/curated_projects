#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os


class Utilities:
    @staticmethod
    def ends_with_slash(string):
        return (string + "/", string)[string.endswith("/")]

    @staticmethod
    def get_page(url):
        import requests
        header = "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US) AppleWebKit/525.19 (KHTML, like Gecko) Chrome/1.0.154.53 Safari/525.19"
        return requests.get(url, headers={'User-Agent': header}).content

    @staticmethod
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

    @staticmethod
    def dict2pd_series(dictionary):
        import pandas as pd
        output = pd.Series()
        for key in dictionary:
            output.at[key] = dictionary[key]
        return output

    @staticmethod
    def multi_core_queue(func, queue):
        import multiprocessing
        pool = multiprocessing.Pool()
        output = pool.map(func, queue)
        pool.close()
        pool.join()
        return output

    @staticmethod
    def single_core_queue(func, queue):
        return [func(i) for i in queue]

    @staticmethod
    def filename_only(string):
        return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))

    @staticmethod
    def get_script_dir():
        import sys
        return Utilities.ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0])))

    @staticmethod
    def compose_sampledatas_dict(dir_name: str):
        import re
        dir_name = Utilities.ends_with_slash(dir_name)
        files_list = os.listdir(dir_name)
        output_dict = {}
        for file_name in files_list:
            if any([file_name.endswith(i) for i in ["csfasta", "fasta", "fa", "fastq", "fq", "gz"]]):
                sample_name = file_name.split("_")[0].strip()
                sample_files = [dir_name + i for i in files_list if sample_name in i]
                sample_files.sort(key=len, reverse=True)
                sample_name = re.sub("[^A-Za-z0-9]+", "_", sample_name)
                sample_name = re.sub("_+", "_", sample_name)
                output_dict[sample_name] = sample_files
        output_dict = dict(sorted(output_dict.items()))
        return output_dict
