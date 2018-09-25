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
    def load_string(file: str):
        with open(file=file, mode="r", encoding="utf-8") as f:
            s = f.read()
        return s

    @staticmethod
    def split_lines(string: str):
        import re
        return Utilities.remove_empty_values([i.strip() for i in re.sub("[\r\n]+", "\n", string).split("\n")])

    @staticmethod
    def load_list(file: str):
        return Utilities.split_lines(Utilities.load_string(file))

    @staticmethod
    def string_to_2d_array(string: str):
        return Utilities.remove_empty_values([[j.strip() for j in i.split("\t")] for i in Utilities.split_lines(string)])

    @staticmethod
    def _2d_array_to_dicts_list(arr: list, names: list):
        if len(arr[0]) != len(names):
            raise ValueError("Cannot parse dictionary: keys number is not equal!")
        out = []
        for row_list in arr:
            counter = 0
            while counter < len(row_list):
                out.append({names[counter]: row_list[counter]})
                counter += 1
        return [[j for j in i] for i in arr]

    @staticmethod
    def load_2d_array(file: str):
        return Utilities.string_to_2d_array(Utilities.load_string(file))

    @staticmethod
    def load_dicts_list(file: str, names: list):
        arr = Utilities.load_2d_array(file)
        if len(arr[0]) != len(names):
            raise ValueError("Cannot parse dictionary: keys number is not equal!")

    @staticmethod
    def dump_string(string: str, file: str):
        with open(file=file, mode="w", encoding="utf-8") as f:
            f.write(string)

    @staticmethod
    def dump_list(lst: list, file: str):
        Utilities.dump_string(string="\n".join([str(i) for i in lst]) + "\n", file=file)

    @staticmethod
    def dump_2d_array(array: list, file: str):
        Utilities.dump_list(lst=["\t".join([str(j) for j in i]) for i in array], file=file)

    @staticmethod
    def _compose_sampledatas_dict(dir_name: str):
        import re
        dir_name = Utilities.ends_with_slash(dir_name)
        files_list = os.listdir(dir_name)
        output_dict = {}
        for file_name in files_list:
            if any([file_name.endswith(i) for i in ["csfasta", "fasta", "fa", "fastq", "fq", "gz"]]):
                sample_name = file_name.split("_")[0].strip()
                sample_extension = sample_name.split(".")[-1]
                sample_files = [dir_name + i for i in files_list if len(re.findall("^{}".format(sample_name), i)) > 0 and i.endswith(sample_extension)][:2]
                sample_files.sort(key=len, reverse=True)
                sample_name = re.sub("[^A-Za-z0-9]+", "_", sample_name)
                sample_name = re.sub("_+", "_", sample_name)
                output_dict[sample_name] = sample_files
        output_dict = dict(sorted(output_dict.items()))
        return output_dict

    @staticmethod
    def create_sampledata(dirs: list, file: str):
        output_dict = {}
        for dir_name in dirs:
            output_dict.update(Utilities._compose_sampledatas_dict(dir_name))
        Utilities.dump_2d_array(array=[[k] + output_dict[k] for k in output_dict], file=file)
