#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

import subprocess
import numpy as np
import re
import os
from auhrbach.klebsiella_infants.ProjectDescriber import ProjectDescriber


# Get the raw reads files
reads_files_dir = "/data1/bio/190405_M01969_0041_000000000-C6B66/Conversion_shotgun/Klebsiella"
# Sort them alphabetically
reads_files_list = np.array(list(sorted(
    [j for j in [i.strip() for i in subprocess.getoutput("ls -d {}/*".format(reads_files_dir)).split("\n")] if
     len(j) > 0])), dtype=str)
# Chop them to chunks with length of 2
reads_files_pairs_2d_array = reads_files_list.reshape(-1, 2)
# Are reads files corresponding to each other? Is their count even?
assert all([i[0] == i[1].replace("_R2_", "_R1_") and i[0].replace("_R1_", "_R2_") == i[1] for i in
            reads_files_pairs_2d_array]) and len(reads_files_list) % 2 == 0
# Get the sample names from reads file names
sampledata_dict = {re.findall("(.+)_S[0-9]{2}_R[1|2]_001.fastq.gz", i[0].split("/")[-1])[0]: list(i) for i in
                   reads_files_pairs_2d_array}
# Prepare sample data for export
sampledata_string = "".join(["{}\t{}\n".format(k, "\t".join(sampledata_dict[k])) for k in sampledata_dict])
# Export sampledata
os.makedirs(ProjectDescriber.directory, exist_ok=True)
raw_sampledata_file = os.path.join(ProjectDescriber.directory, "raw.sampledata")
with open(raw_sampledata_file, mode="w", encoding="utf-8") as file:
    file.write(sampledata_string)

print(raw_sampledata_file)  # /data1/bio/projects/auhrbach/klebsiella_infants/raw.sampledata
