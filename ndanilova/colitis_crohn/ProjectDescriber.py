#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.templates.ProjectDescriberTemplate import ProjectDescriberTemplate


class ProjectDescriber(ProjectDescriberTemplate):
    owner = "ndanilova"
    name = "colitis_crohn"
    directory = "/data1/bio/projects/ndanilova/colitis_crohn/"
    groupdata = ("/data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/all_groups.groupdata",
                 "/data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/remission_vs_escalation.groupdata",
                 "/data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/crohn_vs_colitis.groupdata")
    sampledata = "/data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata"
    mask = "no_hg19"


class DanilovaAwesomeGroupsDefiner:
    def __init__(self):
        # See "SCFAs_from_KEGG" project for details
        self._groupdata_file = "/data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.groupdata"
    def digest_groupdata(self):
        import pandas as pd
        import os
        import subprocess
        main_groupdata_df = pd.read_table(self._groupdata_file, sep="\t", header=0)
        main_groupdata_df_raw_col_names = list(main_groupdata_df)
        main_groupdata_df["all_groups"] = main_groupdata_df["group_id"].apply(self._define_control)
        main_groupdata_df["remission_vs_escalation"] = main_groupdata_df["all_groups"].apply(self._define_remission_vs_escalation)
        main_groupdata_df["crohn_vs_colitis"] = main_groupdata_df["all_groups"].apply(self._define_crohn_vs_colitis)
        output_dir = "{}/group_data_digest/".format(os.path.dirname(self._groupdata_file))
        subprocess.getoutput("rm -rf {}".format(output_dir))
        os.makedirs(output_dir)
        processed_files = []
        for group_digestion_method in [i for i in list(main_groupdata_df) if i not in main_groupdata_df_raw_col_names]:
            digested_groupdata_df = main_groupdata_df.loc[:, ["sample_name", group_digestion_method]]
            output_file = "{}{}.groupdata".format(output_dir, group_digestion_method)
            digested_groupdata_df.to_csv(output_file, sep="\t", header=False, index=False)
            processed_files.append(output_file)
        print("Exported group data digest: \n\"{}\"".format("\", \n\"".join(processed_files)))
    @staticmethod
    def _define_control(s: str):
        s = s.strip()
        if s == "srr":
            return "control"
        return s
    @staticmethod
    def _define_remission_vs_escalation(s: str):
        if s.endswith("esc"):
            return "escalation"
        elif s.endswith("rem"):
            return "remission"
        else:
            return s
    @staticmethod
    def _define_crohn_vs_colitis(s: str):
        return s.split("_")[0].strip()


if __name__ == '__main__':
    describer = ProjectDescriber()
    groupdataDefiner = DanilovaAwesomeGroupsDefiner()
    groupdataDefiner.digest_groupdata()
