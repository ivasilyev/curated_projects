#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os


class Utilities:
    @staticmethod
    def ends_with_slash(string):
        return (string + "/", string)[string.endswith("/")]

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
    def filename_only(string):
        return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))

    @staticmethod
    def get_script_dir():
        import sys
        return Utilities.ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0])))

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
    def dump_list(lst: list, file: str):
        Utilities.dump_string(string="\n".join([str(i) for i in lst]) + "\n", file=file)

    @staticmethod
    def dump_2d_array(array: list, file: str):
        Utilities.dump_list(lst=["\t".join([str(j) for j in i]) for i in array], file=file)

    @staticmethod
    def flatten_2d_array(array: list):
        return [j for i in array for j in i]

    @staticmethod
    def safe_findall(pattern, string, idx: int = 0):
        import re
        try:
            return re.findall(pattern, string)[idx]
        except IndexError:
            print("Warning! Can't find the regex pattern '{}' within the string: '{}'".format(pattern, string))
            return ""

    # File processing methods

    @staticmethod
    def load_string(file: str):
        with open(file=file, mode="r", encoding="utf-8") as f:
            s = f.read()
        return s

    @staticmethod
    def dump_string(string: str, file: str):
        with open(file=file, mode="w", encoding="utf-8") as f:
            f.write(string)

    @staticmethod
    def scan_whole_dir(dir_name: str):
        out = []
        for root, dirs, files in os.walk(dir_name):
            for file in files:
                out.append(os.path.join(root, file))
        return out

    @staticmethod
    def concatenate_files(*source_files, target_file):
        if len(source_files) < 2:
            raise ValueError("Not enough files (at least 2 needed): '{}'".format(source_files))
        with open(target_file, mode="w", encoding="utf-8") as target:
            for source_file in source_files:
                with open(source_file, mode="r", encoding="utf-8") as source:
                    for line in source:
                        target.write(line)
                    target.write("\n")
                    source.close()
            target.close()

    # Pandas methods

    @staticmethod
    def dict2pd_series(dictionary):
        import pandas as pd
        output = pd.Series()
        for key in dictionary:
            output.at[key] = dictionary[key]
        return output

    @staticmethod
    def merge_pd_series_list(series: list):
        import pandas as pd
        return pd.concat(series, axis=1, sort=False).transpose()

    @staticmethod
    def left_merge(df0, df1, *args):
        import pandas as pd
        assert isinstance(df0, pd.DataFrame) and isinstance(df1, pd.DataFrame)
        df1_uniq = [i for i in list(df1) if i not in list(df0)]
        return df0.merge(df1.loc[:, list(args) + df1_uniq], on=args, how="left")

    @staticmethod
    def combine_duplicate_rows(df, index_col_name):
        def _combine(*args):
            strings = set([j for j in [str(i).strip() for i in args] if len(j) > 0])
            return ";".join(strings)
        import pandas as pd
        for dupe_id in df[df[index_col_name].duplicated()][index_col_name].values.tolist():
            subdf = df.loc[df[index_col_name] == dupe_id, :].copy()
            df = df.loc[df[index_col_name] != dupe_id, :]
            accumulator = pd.Series()
            for idx in range(0, len(subdf.index.values.tolist())):
                row = subdf.iloc[idx].fillna("")
                if len(accumulator) == 0:
                    accumulator = row
                else:
                    accumulator = accumulator.combine(row, _combine)
            df = pd.concat([df, pd.DataFrame(accumulator).transpose()], axis=0, ignore_index=True)
        df.sort_values(index_col_name, inplace=True)
        return df

    # Queue processing methods

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

    # Web methods

    @staticmethod
    def get_page(url):
        import requests
        header = "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US) AppleWebKit/525.19 (KHTML, like Gecko) Chrome/1.0.154.53 Safari/525.19"
        return requests.get(url, headers={'User-Agent': header}).content

    @staticmethod
    def scrap_links_from_web_page(url: str) -> list:
        import bs4
        import lxml
        soup = bs4.BeautifulSoup(Utilities.get_page(url), "lxml")
        root = "/".join(url.split("/")[:-1])
        # Scrap the web page roughly and get all required links
        raw_links = [a["href"].strip() for a in soup.find_all("a", href=True)]
        return [i if "://" in i else "{}/{}".format(root, i) for i in raw_links]

    @staticmethod
    def download_file(url, out_dir):
        import subprocess
        from time import sleep
        _RETRIES_LEFT = 5
        _SLEEP_SECONDS = 3
        url = url.strip()
        out_dir = os.path.normpath(out_dir.strip())
        assert len(url) > 0 and len(out_dir) > 0
        os.makedirs(out_dir, exist_ok=True)
        while _RETRIES_LEFT > 0:
            out_file = os.path.join(out_dir, url.split("/")[-1])
            print(subprocess.getoutput("curl -fsSL {} -o {}".format(url, out_file)))
            sleep(_SLEEP_SECONDS)
            if os.path.isfile(out_file):
                print("Download finished: '{}'".format(out_file))
                return out_file
            _RETRIES_LEFT -= 1
            print("Warning! Failed download: '{}'. Retries left: {}".format(url, _RETRIES_LEFT))
        return ""
