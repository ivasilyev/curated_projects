#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -e DISPLAY=$DISPLAY -it ${DOCKER_IMAGE_NAME} python3
"""

import pandas as pd


class PivotTableAnnotator:
    def __init__(self, pivot_file, annotation_file):
        pivot_df = pd.read_table(pivot_file, sep='\t', header=0, engine='python')
        annotation_df = pd.read_table(annotation_file, sep='\t', header=0, engine='python')
        self.annotated_df = pd.merge(annotation_df, pivot_df, on="reference_id", how='outer')
    def export_annotated_pivot(self, output_file):
        self.annotated_df.to_csv(path_or_buf=output_file, sep='\t', header=True, index=False)
    @staticmethod
    def _draw_boxplot(left_list, right_list, left_name, right_name, title, output_file):
        import matplotlib as mpl
        ## agg backend is used to create plot as a .png file
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        data_to_plot = [left_list, right_list]
        # Create a figure instance
        fig = plt.figure(1, figsize=(9, 6))
        # Create an axes instance
        ax = fig.add_subplot(111)
        # Create the boxplot
        def draw_plot(data, edge_color, fill_color):
            bp = ax.boxplot(data, patch_artist=True)
            for element in ['boxes', 'whiskers', 'fliers', 'means', 'caps']:
                plt.setp(bp[element], color=edge_color)
            plt.setp(bp['medians'], color='black')
            for patch in bp['boxes']:
                patch.set(facecolor=fill_color)
        draw_plot(data_to_plot, 'red', 'tan')
        plt.title(title, y=1.08)
        ax.set_xticklabels([left_name, right_name], rotation=0, fontsize=8)
        ax.set_xlabel('Groups')
        ax.set_ylabel('Values')
        # Save the figure
        fig.savefig(output_file, bbox_inches='tight')
        # Closes the current figure
        plt.close()
    def export_group_reports(self, output_dir, multi_test="fdr_bh", single_test="u-test", base_col_name="reference_id", annotation_col_name="former_id"):
        import os
        output_dir = (output_dir + "/", output_dir)[output_dir.endswith("/")]
        bool_suffix = "_is_rejected_by_{m}_for_{s}".format(m=multi_test, s=single_test)
        bool_col_names_list = [i for i in list(self.annotated_df) if i.endswith(bool_suffix)]
        for bool_col_name in bool_col_names_list:
            comparison_string = bool_col_name.replace(bool_suffix, "")
            comparison_list = comparison_string.split("_vs_")
            left_group_name = comparison_list[0]
            right_group_name = comparison_list[1]
            comparison_df = self.annotated_df.loc[self.annotated_df[bool_col_name] == True].set_index(base_col_name)
            if len(comparison_df) > 0:
                left_group_pivot_df = comparison_df.loc[:, [i for i in comparison_df if i.split("_")[0] == left_group_name and i.split("_")[1].startswith("/")]]
                right_group_pivot_df = comparison_df.loc[:, [i for i in comparison_df if i.split("_")[0] == right_group_name and i.split("_")[1].startswith("/")]]
                comparison_dir = "{a}{b}/".format(a=output_dir, b=comparison_string)
                os.makedirs(comparison_dir, exist_ok=True)
                comparison_df.to_csv(path_or_buf="{}pivot.tsv".format(comparison_dir), sep='\t', header=True, index=True)
                left_group_pivot_df.to_csv(path_or_buf="{a}raw_{b}.tsv".format(a=comparison_dir, b=left_group_name), sep='\t', header=True, index=True)
                right_group_pivot_df.to_csv(path_or_buf="{a}raw_{b}.tsv".format(a=comparison_dir, b=right_group_name), sep='\t', header=True, index=True)
                boxplots_dir = "{}boxplots/".format(comparison_dir)
                os.makedirs(boxplots_dir, exist_ok=True)
                for base_id in comparison_df.index.tolist():
                    annotated_id = comparison_df.loc[base_id, annotation_col_name]
                    # NaN is float type
                    if not isinstance(annotated_id, str):
                        annotated_id = comparison_df.loc[base_id, "former_id"]
                    self._draw_boxplot(left_list=left_group_pivot_df.loc[base_id].values.tolist(),
                                       right_list=right_group_pivot_df.loc[base_id].values.tolist(),
                                       left_name=left_group_name,
                                       right_name=right_group_name,
                                       title=annotated_id,
                                       output_file="{a}{b}.png".format(a=boxplots_dir, b=base_id))
