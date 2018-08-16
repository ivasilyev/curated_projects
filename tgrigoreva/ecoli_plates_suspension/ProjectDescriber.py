#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""


class ProjectDescriber:
    owner = "tgrigoreva"
    name = "ecoli_plates_suspension"
    directory = "/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/"
    groupdata = "/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/raw.groupdata"
    sampledata = "/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/raw.sampledata"
    mask = "no_hg19"


class ProjectEvaluator:
    @staticmethod
    def evaluate_sampledata():
        import os
        import subprocess
        import pandas as pd
        from meta.scripts.utilities import filename_only
        import re
        #
        df = pd.DataFrame(columns=["sample_name", "sample_path"])
        for dir_mask in ["/data2/bio/ecoli_komfi/raw_reads/*", "/data2/bio/ecoli_komfi/raw_reads2/*"]:
            data_1 = [i.strip() for i in subprocess.getoutput("ls -d {}R1*.fastq* | sort".format(dir_mask)).split("\n")]
            data_12 = ["{a}\t{b}".format(a=i, b=i.replace("R1", "R2")) if os.path.isfile(i.replace("R1", "R2")) else "" for i in data_1]
            sample_names_list = [re.sub("_S.*$", "", filename_only(i)) for i in data_1]
            df = pd.concat([df, pd.DataFrame.from_dict({"sample_name": sample_names_list, "sample_path": data_12})], axis=0, ignore_index=True)
        #
        df["group_id"] = "group_id"
        #
        os.makedirs(ProjectDescriber.directory, exist_ok=True)
        df.loc[:, ["sample_name", "group_id"]].to_csv(ProjectDescriber.groupdata, sep='\t', index=False, header=False)
        df.loc[:, ["sample_name", "sample_path"]].to_csv(ProjectDescriber.sampledata, sep='\t', index=False, header=False)
        subprocess.getoutput("sed -i 's|\"||g' {}".format(ProjectDescriber.sampledata))  # Tab-containing columns items are flanked by '"'
