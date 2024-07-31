#%%

import os
from meta.utils.file_system import find_by_regex
from meta.utils.file_system import filename_only
from meta.utils.pandas import extract_tables_from_zip_into_excel

pipeline_output_dir = "/data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/qiime2-picrust2-pipeline-out"
qzv_files = find_by_regex(".*\.qzv$", pipeline_output_dir)
excel_files_dir = os.path.join(pipeline_output_dir, "xlsx")

#%%

for qzv_file in qzv_files:
    qzv_filename = filename_only(qzv_file)
    excel_file = os.path.join(excel_files_dir, f"{qzv_filename}.xlsx")
    extract_tables_from_zip_into_excel(qzv_file, excel_file)
