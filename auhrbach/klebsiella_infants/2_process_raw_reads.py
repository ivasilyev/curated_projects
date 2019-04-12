# -*- coding: utf-8 -*-

"""
# From another console:
export IMG=ivasilyev/spades_cutadapt:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

import subprocess
import re
import os


def multi_core_queue(func, queue):
    import multiprocessing
    pool = multiprocessing.Pool()
    output = pool.map(func, queue)
    pool.close()
    pool.join()
    return output


def sampledata_line_to_cutadapt(s: str):
    if len(s.strip()) == 0:
        return
    try:
        sample_name, raw_file_1, raw_file_2 = [i.strip() for i in s.split("\t")]
    except ValueError:
        raise ValueError("Cannot process the sample data line: '{}'".format(s))
    cutadapt_dir = os.path.join(os.path.dirname(raw_file_1), "cutadapt")
    os.makedirs(cutadapt_dir, exist_ok=True)
    _ADAPTER = "AGATCGGAAGAG"
    trimmed_file_1, trimmed_file_2 = [
        os.path.join(cutadapt_dir, "{}_cutadapt.{}.fastq".format(sample_name, re.findall("_R([\d])_", i)[0])) for i
        in [raw_file_1, raw_file_2]]
    log_file = os.path.join(cutadapt_dir, "{}_cutadapt.log".format(sample_name))
    for _f in [trimmed_file_1, trimmed_file_2, log_file]:
        if os.path.exists(_f):
            os.remove(_f)
    cmd = "cutadapt -a {ad} -A {ad} -m 50 -o {t1} -p {t2} {r1} {r2}".format(ad=_ADAPTER, t1=trimmed_file_1,
                                                                            t2=trimmed_file_2, r1=raw_file_1,
                                                                            r2=raw_file_2)
    try:
        log = subprocess.getoutput(cmd)
    except PermissionError:
        raise ValueError("Permission denied, please run `sudo chmod -R 777 {}`".format(os.path.dirname(raw_file_1)))
    with open(log_file, mode="w", encoding="utf-8") as f:
        f.write(log)
    return {"sample_name": sample_name, "trimmed_file_1": trimmed_file_1, "trimmed_file_2": trimmed_file_2}


def single_core_queue(func, queue):
    output = [func(i) for i in queue]
    return output


def sampledata_dict_to_spades(d: dict):
    # Check all keys and values exist
    assert all([d.get(i) for i in ["sample_name", "trimmed_file_1", "trimmed_file_2"]])
    spades_dir = os.path.join(os.path.dirname(os.path.dirname(d.get("trimmed_file_1"))), "spades", d.get("sample_name"))
    subprocess.getoutput("rm -rf {}".format(spades_dir))
    os.makedirs(spades_dir)
    log_file = os.path.join(spades_dir, "{}_spades.log".format(d.get("sample_name")))
    cmd = "python3 /opt/spades/bin/spades.py --careful -o {out} -1 {t1} -2 {t2}".format(out=spades_dir,
                                                                                        t1=d.get("trimmed_file_1"),
                                                                                        t2=d.get("trimmed_file_2"),
                                                                                        log=log_file)
    log = subprocess.getoutput(cmd)
    with open(log_file, mode="w", encoding="utf-8") as f:
        f.write(log)
    return os.path.join(spades_dir, "contigs.fasta")


raw_sampledata_file = "/data1/bio/projects/auhrbach/klebsiella_infants/raw.sampledata"
with open(raw_sampledata_file, mode="r", encoding="utf-8") as file:
    raw_sampledata_list = [i.strip() for i in file.read().split("\n") if len(i.strip()) > 0]

trimmed_sampledatas_list = multi_core_queue(sampledata_line_to_cutadapt, raw_sampledata_list)
trimmed_sampledata_file = os.path.join(os.path.dirname(raw_sampledata_file), "trimmed.sampledata")
with open(trimmed_sampledata_file, mode="w", encoding="utf-8") as file:
    file.write("".join(
        sorted(["{}\t{}\t{}\n".format(i.get("sample_name"), i.get("trimmed_file_1"), i.get("trimmed_file_2")) for i in
                trimmed_sampledatas_list])))

assembled_genomes_list = single_core_queue(sampledata_dict_to_spades, trimmed_sampledatas_list)
