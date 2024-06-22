#%%

import os
from collections import deque
from meta.utils.file_system import find_by_regex, copy, to_7z


pipeline_output_dir = "/data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/qiime2-picrust2-pipeline"
target_dir = os.path.join(os.path.dirname(pipeline_output_dir), "qiime2-picrust2-pipeline-out")

copying_dicts = deque([
    {
        "source": i,
        "destination": os.path.join(target_dir, "dada2", "qzv", f"dada2_{os.path.basename(i)}"),
    }
    for i in find_by_regex(".+qiime2.+dada2.+\.qzv$", pipeline_output_dir)
])

#%%

copying_dicts.extend([
    {
        "source": i,
        "destination": os.path.join(target_dir, "vsearch", "qzv", f"vsearch_{os.path.basename(i)}"),
    }
    for i in find_by_regex(".+qiime2.+vsearch.+\.qzv$", pipeline_output_dir)
])
copying_dicts

#%%

for copying_dict in copying_dicts:
    copy(**copying_dict)

#%%

compressing_dicts = deque()
asv_file = find_by_regex(".+qiime2.+dada2.+ASV_with_taxa.tsv$", pipeline_output_dir)[0]
compressing_dicts.append({
    "source": asv_file,
    "destination": os.path.join(target_dir, "dada2", "tsv", f"dada2_{os.path.basename(asv_file)}.7z"),
})
compressing_dicts

#%%

otu_normalized_file = find_by_regex(".+qiime2.+vsearch.+OTU_normalized_with_taxa.tsv$", pipeline_output_dir)[0]
compressing_dicts.append({
    "source": otu_normalized_file,
    "destination": os.path.join(target_dir, "vsearch", "tsv", f"vsearch_{os.path.basename(otu_normalized_file)}.7z"),
})
compressing_dicts

#%%

for picrust2_output_filename in (
    "EC_pred_metagenome_unstrat_described",
    "KO_pred_metagenome_unstrat_described",
    "KO_pred_metagenome_unstrat_described_annotated",
    "pathways_abun_unstrat_described",
):
    print(f"Add '{picrust2_output_filename}'")
    picrust2_output_file = find_by_regex(f".+picrust2.+{picrust2_output_filename}\.tsv$", pipeline_output_dir)[0]
    compressing_dicts.append({
        "source": picrust2_output_file,
        "destination": os.path.join(target_dir, "picrust2", "tsv", f"{os.path.basename(picrust2_output_file)}.7z")
    })

compressing_dicts

#%%

for compressing_dict in compressing_dicts:
    to_7z(**compressing_dict)
