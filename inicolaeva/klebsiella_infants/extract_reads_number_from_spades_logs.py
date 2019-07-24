#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull
python3
"""

import os
import re
import shutil
import subprocess
import pandas as pd
from meta.scripts.Utilities import Utilities
from inicolaeva.klebsiella_infants.ProjectDescriber import ProjectDescriber


log_dir = "/data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/log/2019-06-28-11-07-17"
spades_log_files = [i for i in Utilities.scan_whole_dir(log_dir) if os.path.basename(i).startswith("spades")]
spades_log_dict = [{"sample_name": Utilities.safe_findall("spades_(.+).log", i), "spades_log_file": i}
                   for i in spades_log_files]
i = '/data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/log/2019-06-28-11-07-17/spades_Kleb102.log'

[j for j in Utilities.load_list(i) if all(m in j for m in ["Total", "reads processed"])]






