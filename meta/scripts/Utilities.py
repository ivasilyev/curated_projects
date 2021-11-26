#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os


class Utilities:
    # File system based methods
    @staticmethod
    def is_file_valid(file: str, report: bool = True):
        if not os.path.exists(file):
            if report:
                print("Not found: '{}'".format(file))
            return False
        if not os.path.isfile(file):
            if report:
                print("Not a file: '{}'".format(file))
            return False
        if os.path.getsize(file) == 0:
            if report:
                print("Empty file: '{}'".format(file))
            return False
        return True

    @staticmethod
    def scan_whole_dir(dir_name: str):
        out = []
        for root, dirs, files in os.walk(dir_name):
            for file in files:
                out.append(os.path.join(root, file))
        return sorted(out)

    @staticmethod
    def find_file_by_tail(dir_name: str, tail: str, multiple: bool = False):
        files = [i for i in Utilities.scan_whole_dir(dir_name) if i.endswith(tail)]
        if len(files) == 0:
            return ""
        if multiple:
            return files
        return files[0]

    @staticmethod
    def filename_only(s: str):
        return os.path.splitext(os.path.basename(s))[0]

    @staticmethod
    def ends_with_slash(string):
        return (string + "/", string)[string.endswith("/")]

    @staticmethod
    def get_script_dir():
        import sys
        return Utilities.ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0])))

    # I/O methods

    @staticmethod
    def load_string(file: str):
        with open(file=file, mode="r", encoding="utf-8") as f:
            s = f.read()
            f.close()
        return s

    @staticmethod
    def dump_string(string: str, file: str):
        os.makedirs(os.path.dirname(file), exist_ok=True)
        with open(file=file, mode="w", encoding="utf-8") as f:
            f.write(string)
            f.close()

    @staticmethod
    def load_dict(file: str):
        import json
        return json.loads(Utilities.load_string(file))

    @staticmethod
    def dump_dict(d: dict, file: str, **kwargs):
        _kwargs = dict(indent=4, sort_keys=False)
        if len(kwargs.keys()) > 0:
            _kwargs.update(kwargs)
        import json
        return Utilities.dump_string(json.dumps(d, **_kwargs), file)

    @staticmethod
    def load_list(file: str):
        return Utilities.split_lines(Utilities.load_string(file))

    @staticmethod
    def dump_list(lst: list, file: str):
        Utilities.dump_string(string="\n".join([str(i) for i in lst]) + "\n", file=file)

    @staticmethod
    def load_2d_array(file: str):
        return Utilities.string_to_2d_array(Utilities.load_string(file))

    @staticmethod
    def dump_2d_array(array: list, file: str):
        Utilities.dump_list(lst=["\t".join([str(j) for j in i]) for i in array], file=file)

    @staticmethod
    def load_dicts_list(file: str, names: list):
        arr = Utilities.load_2d_array(file)
        if len(arr[0]) != len(names):
            raise ValueError("Cannot parse dictionary: keys number is not equal!")

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

    @staticmethod
    def decompress_file(file_name: str, remove_file: bool = True):
        import subprocess
        _ARCHIVE_EXTENSIONS = ("tar", "gz", "bz2", "zip")
        cmd = ""
        if not any(file_name.endswith(".{}".format(i)) for i in _ARCHIVE_EXTENSIONS):
            print("Nothing to extract: '{}'".format(file_name))
        if file_name.endswith(".tar.gz"):
            cmd = "tar -xvzf {i} -C {o}/"
        elif file_name.endswith(".tar.bz2"):
            cmd = "tar -xvjf {i} -C {o}/"
        elif file_name.endswith(".tar"):
            cmd = "tar -xvf {i} -C {o}/"
        elif file_name.endswith(".gz"):
            cmd = "cd {o} && gzip -d {i}"
        elif file_name.endswith(".zip"):
            cmd = "unzip {i} -d {o}/"
        elif file_name.endswith(".rar"):
            cmd = "unrar e {i} {o}/"
        print(subprocess.getoutput(cmd.format(i=file_name, o=os.path.normpath(os.path.dirname(file_name)))))
        print("Extracting completed: '{}'".format(file_name))
        if os.path.isfile(file_name) and remove_file:
            print("Removing file: '{}'".format(file_name))
            os.remove(file_name)

    # System methods

    @staticmethod
    def count_elapsed_seconds(t):
        from time import perf_counter
        return f"{perf_counter() - t :.3f} s."

    @staticmethod
    def get_time():
        from datetime import datetime
        now = datetime.now()
        output_list = []
        for time_unit in [now.year, now.month, now.day, now.hour, now.minute, now.second]:
            time_unit = str(time_unit)
            if len(time_unit) < 2:
                time_unit = '0' + time_unit
            output_list.append(time_unit)
        return '-'.join(output_list)

    # Primitive processing methods

    @staticmethod
    def safe_findall(pattern, string, idx: int = 0, report: bool = False):
        import re
        try:
            return re.findall(pattern, string)[idx]
        except IndexError:
            if report:
                print("Warning! Can't find the regex pattern '{}' within the string: '{}'".format(
                    pattern, string))
            return ""

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
    def split_lines(string: str):
        import re
        out = [i.strip() for i in re.sub("[\r\n]+", "\n", string).split("\n")]
        return Utilities.remove_empty_values(out)

    @staticmethod
    def string_to_2d_array(string: str):
        out = [[j.strip() for j in i.split("\t")] for i in Utilities.split_lines(string)]
        return Utilities.remove_empty_values(out)

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
    def flatten_2d_array(array: list):
        return [j for i in array for j in i]

    @staticmethod
    def split_list_by_chunk_length(list_: list, n: int):
        return [list_[i:i + n] for i in range(0, len(list_), n)]

    @staticmethod
    def split_list_by_chunk_number(list_: list, n: int):
        from math import ceil
        return Utilities.split_list_by_chunk_length(list_, ceil(len(list_) / n))

    @staticmethod
    def get_most_similar_word_from_list(word: str, words: list):
        from difflib import SequenceMatcher
        _words = [i for i in sorted(set(words)) if i != word]
        sorted_ratios = sorted(
            {i: SequenceMatcher(a=word, b=i).ratio() for i in _words}.items(),
            key=lambda x: x[1],
            reverse=True
        )
        return sorted_ratios[0][0]

    @staticmethod
    def get_most_similar_word_pairs(words: list):
        _words = sorted(set(words))
        out = []
        for word in _words:
            word_pair = sorted([word, Utilities.get_most_similar_word_from_list(word, _words)])
            if word_pair not in out:
                out.append(word_pair)
        return out

    @staticmethod
    def filtered_product(lists: list):
        from itertools import product
        out = []
        for sub_list in list(product(*lists)):
            sub_out = []
            is_sub_good = True
            for item in sub_list:
                is_sub_good = item not in sub_out
                if is_sub_good:
                    sub_out.append(item)
                else:
                    break
            if is_sub_good and sorted(sub_out) not in [sorted(i) for i in out]:
                out.append(tuple(sub_out))
        return out

    @staticmethod
    def is_sequence_file(s: str):
        return any(s.endswith(i) for i in [
            ".fasta", ".fastq.gz", ".fastq", ".fq.gz", "fq"
    ])

    # Biopython methods

    @staticmethod
    def parse_sequences(file: str, type_: str = "fasta"):
        from Bio import SeqIO
        if type_ == "fastq_gz":
            import gzip
            with gzip.open(file, "rt") as f:
                records = list(SeqIO.parse(f, "fastq"))
                f.close()
        else:
            with open(file, mode="r", encoding="utf-8") as f:
                records = list(SeqIO.parse(f, type_))
                f.close()
        out = Utilities.remove_empty_values(sorted(records, key=lambda x: len(x), reverse=True))
        return out

    @staticmethod
    def remove_duplicate_sequences(records: list):
        out = []
        sequences = []
        for record in records:
            if record.seq not in sequences:
                sequences.append(record.seq)
                out.append(record)
        return out

    @staticmethod
    def randomize_gene_slice(record, size: int = 20000):
        """
        :param record: SeqRecord

        :param size: int
        The typical gene is about 1000 bp in length:
        http://bioscience.jbpub.com/cells/MBIO137.aspx
        The default slicing will return a chunk containing ~20 genes.

        :return: trimmed SeqRecord
        """
        from copy import deepcopy
        from random import randint
        gene_length = len(record)
        if gene_length <= size:
            return record
        start = randint(0, gene_length - size)
        end = start + size
        record_ = deepcopy(record)
        record_.seq = record_.seq[start:end]
        return record_

    @staticmethod
    def count_reads_statistics(reads_file: str, type_: str = "fasta", prefix: str = "",
                               kwargs: dict = None) -> dict:
        import statistics
        from Bio.SeqUtils import GC
        seq_records = Utilities.parse_sequences(reads_file, type_)
        total_sequence = "".join([str(i.seq) for i in seq_records])
        out = {
            "reads_file": reads_file,
            "reads_number": len(seq_records),
            "largest_read_bp": len(seq_records[0]),
            "smallest_read_bp": len(seq_records[-1]),
            "total_reads_bp": len(total_sequence),
            "mean_reads_bp": statistics.mean([len(i) for i in seq_records]),
            "median_reads_bp": statistics.median([len(i) for i in seq_records]),
            "gc_percentage": GC(total_sequence),
        }
        if isinstance(kwargs, dict) and len(kwargs.keys()) > 0:
            out.update(kwargs)
        if len(prefix) > 0:
            return {"{}_{}".format(prefix, k): out[k] for k in out.keys()}
        return out

    @staticmethod
    def count_assembly_statistics(assembly_file: str, type_: str = "fasta",
                                  prefix: str = "", kwargs: dict = None) -> dict:
        import statistics
        from Bio.SeqUtils import GC
        seq_records = Utilities.parse_sequences(assembly_file, type_=type_)
        total_sequence = "".join([str(i.seq) for i in seq_records])
        out = {
            "assembly_file": assembly_file,
            "contigs_number": len(seq_records),
            "largest_contig_bp": len(seq_records[0]),
            "smallest_contig_bp": len(seq_records[-1]),
            "total_contigs_bp": len(total_sequence),
            "mean_contigs_bp": statistics.mean([len(i) for i in seq_records]),
            "median_contigs_bp": statistics.median([len(i) for i in seq_records]),
            "gc_percentage": GC(total_sequence),
        }
        contigs_records = []
        bp50 = round(len(total_sequence) * 0.5)
        bp90 = round(len(total_sequence) * 0.9)
        n50_supplied = False
        n90_supplied = False
        for idx, seq_record in enumerate(seq_records):
            contigs_records.append(seq_record)
            if sum([len(i) for i in contigs_records]) > bp50 and not n50_supplied:
                out["l50"] = idx + 1
                out["n50"] = len(seq_records[idx])
                out["d50"] = len(seq_records[idx + 1])
                n50_supplied = True
            if sum([len(i) for i in contigs_records]) > bp90 and not n90_supplied:
                out["l90"] = idx + 1
                out["n90"] = len(seq_records[idx])
                n90_supplied = True
        if isinstance(kwargs, dict) and len(kwargs.keys()) > 0:
            out.update(kwargs)
        if len(prefix) > 0:
            return {"{}_{}".format(prefix, k): out[k] for k in out.keys()}
        return out

    @staticmethod
    def count_assembly_coverages(raw_reads_length_sum: int, assembly_length: int,
                                 reference_length: int) -> dict:
        """
        :param raw_reads_length_sum: The number of bases sequenced
        :param assembly_length: The bases that were placed in the final assembly
        :param reference_length: The expected genome size
        :return: dict

        From NCBI template ('Template_GenomeBatch.11700383121d.xlsx'):
        The estimated base coverage across the genome, eg 12x.
        This can be calculated by dividing the number of bases sequenced by the expected genome size
        and multiplying that by the percentage of bases that were placed in the final assembly.
        More simply it is the number of bases sequenced divided by the expected genome size.
        """

        def _process_float(x):
            return "{0:.2f}x".format(x)

        assembled_reads_rate = assembly_length / raw_reads_length_sum
        expected_assembly_coverage = raw_reads_length_sum / reference_length
        real_assembly_coverage = raw_reads_length_sum * assembled_reads_rate / reference_length
        return dict(assembled_reads_rate=assembled_reads_rate,
                    expected_assembly_coverage=_process_float(expected_assembly_coverage),
                    real_assembly_coverage=_process_float(real_assembly_coverage))

    # Pandas methods

    @staticmethod
    def load_tsv(table, col_names: list = None):
        import pandas as pd
        if col_names:
            return pd.read_csv(table, encoding="utf-8", sep="\t", header="infer", names=col_names)
        return pd.read_csv(table, encoding="utf-8", sep="\t", header=0)

    @staticmethod
    def dump_tsv(df, table_file: str, col_names: list = None, reset_index: bool = False):
        import pandas as pd
        assert isinstance(df, pd.DataFrame)
        _df = df.copy()
        os.makedirs(os.path.dirname(table_file), exist_ok=True)
        if col_names is not None and len(col_names) > 0:
            _df = _df.loc[:, col_names]
        if reset_index:
            _df.reset_index(inplace=True)
        _df.to_csv(table_file, encoding="utf-8", sep="\t", index=False, header=True)

    @staticmethod
    def dict2pd_series(dictionary, sort_keys: bool = False):
        import pandas as pd
        output = pd.Series()
        keys = list(dictionary.keys())
        if sort_keys:
            keys = sorted(keys)
        for key in keys:
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

    @staticmethod
    def get_n_majors_from_df(df, col_name: str, n: int = 10, drop_nulls: bool = True):
        """
        :param df: Pandas DataFrame object
        :param col_name: Name of the value column to sort
        :param n: number of rows without to return without the "Other" row
        :param drop_nulls: Only keep values more than 0
        :return: Pandas DataFrame with given number of rows and extra row containing sum of all non-included rows called
                 "Other"
        """
        import pandas as pd
        if not isinstance(df, pd.DataFrame):
            raise ValueError("The first argument must be a Pandas DataFrame object")
        majors_df = df.loc[:, [col_name]].sort_values(col_name, ascending=False).head(n=n)
        if drop_nulls:
            majors_df = majors_df.loc[majors_df[col_name].astype(float) > 0]
        others_df = df.loc[:, [col_name]].sort_values(col_name, ascending=False).tail(n=df.shape[0] - n)
        if others_df.astype(float).sum().values.tolist()[0] == 0.0:
            return majors_df
        out = pd.concat([majors_df, pd.DataFrame(others_df.sum().rename({col_name: "Other"}).rename(col_name))], axis=0)
        out.index.name = df.index.name
        out.columns.name = df.columns.name
        return out

    # Function handling methods

    @staticmethod
    def randomize_sleep(min_: int = 30, max_: int = 120):
        from time import sleep
        from random import randint
        sleep(randint(min_, max_))

    @staticmethod
    def get_caller_name():
        import inspect
        return str(inspect.stack()[1][3])

    @staticmethod
    def attempt_func(func, args):
        _ATTEMPTS = 5
        attempt = 1
        while attempt <= _ATTEMPTS:
            try:
                if any(isinstance(args, i) for i in (list, tuple)):
                    return func(*args)
                if any(isinstance(args, i) for i in (dict,)):
                    return func(**args)
            except Exception as e:
                print("Caught exception for attempt {}: `{}`".format(attempt, e))
                attempt += 1
                Utilities.randomize_sleep()
        print("Exceeded number of attempts for the function: '{}'".format(func.__name__))
        return

    # Queue processing methods

    @staticmethod
    def single_core_queue(func, queue) -> list:
        return [func(i) for i in queue]

    @staticmethod
    def multi_core_queue(func, queue: list, processes: int = 0, async_: bool = False) -> list:
        import multiprocessing as mp
        if processes == 0:
            processes = mp.cpu_count()
        pool = mp.Pool(processes=processes)
        if async_:
            result = pool.map_async(func, queue)
        else:
            result = pool.map(func, queue)
        pool.close()
        pool.join()
        if async_:
            return result.get()
        return result

    @staticmethod
    def wrapper(d: dict):
        """
        An ultimate wrapper
        :param d: Dictionary {'func': a function or a method, 'args': list, 'kwargs': dict}
        :return: The result of evaluating func(*args, **kwargs)
        """
        if "args" not in d.keys():
            d["args"] = []
        if "kwargs" not in d.keys():
            d["kwargs"] = dict()
        return d["func"](*d["args"], **d["kwargs"])

    # Web-based methods

    @staticmethod
    def get_page(url: str, header: str = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                                         "AppleWebKit/537.36 (KHTML, like Gecko) "
                                         "Chrome/90.0.4430.19 Safari/537.36"):
        import requests
        if header is not None:
            return requests.get(url, headers={'User-Agent': header}).content
        return requests.get(url).content

    @staticmethod
    def get_soup(*args, **kwargs):
        import bs4
        import lxml
        return bs4.BeautifulSoup(Utilities.get_page(*args, **kwargs), "lxml")

    @staticmethod
    def scrap_links_from_web_page(url: str) -> list:
        soup = Utilities.get_soup(url)
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
            if os.path.isfile(out_file):
                print("Download finished: '{}'".format(out_file))
                return out_file
            _RETRIES_LEFT -= 1
            sleep(_SLEEP_SECONDS)
            print("Warning! Failed download: '{}'. Retries left: {}".format(url, _RETRIES_LEFT))
        print("Exceeded URL download limits: '{}'".format(url))
        return ""
