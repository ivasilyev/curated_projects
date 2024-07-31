
import os
from argparse import ArgumentParser
from collections import deque
from meta.utils.file_system import filename_only, get_file_extension
from meta.utils.file_system import copy, find_by_regex, to_7z
from meta.utils.pandas import dump_tsv, extract_tables_from_zip_into_excel, load_tsv
from meta.utils.qiime import OTU_COLUMN_NAME


def annotate_qiime2_asv_table(pipeline_output_dir: str):
    asv_file_regex = ".+qiime2.+dada2.+ASV.tsv$"
    mapper_file_regex = ".+qiime2.+dada2.+ASV_confidences.tsv$"
    asv_raw_file = find_by_regex(asv_file_regex, pipeline_output_dir)[0]
    asv_raw_df = load_tsv(asv_raw_file)
    asv_df = asv_raw_df.drop(["taxonomy"], axis=1)
    otu_asv_mapper_file = find_by_regex(mapper_file_regex, pipeline_output_dir)[0]
    otu_asv_mapper_df = load_tsv(otu_asv_mapper_file)
    asv_annotated_df = otu_asv_mapper_df.merge(
        asv_df,
        how="right",
        on=OTU_COLUMN_NAME
    )
    asv_merged_file_basename = f"{filename_only(asv_raw_file)}_with_taxa{get_file_extension(asv_raw_file)}"
    asv_merged_file = os.path.join(os.path.dirname(asv_raw_file), asv_merged_file_basename)
    print(asv_merged_file)
    dump_tsv(asv_annotated_df, asv_merged_file)


def annotate_picrust2_ko_table(pipeline_output_dir, kegg_reference_file: str):
    ko_file_regex = ".+picrust2.+KO_pred_metagenome_unstrat_described.tsv$"
    ko_sample_file = find_by_regex(ko_file_regex, pipeline_output_dir)[0]
    kegg_reference_df = load_tsv(kegg_reference_file)
    ko_sample_df = load_tsv(ko_sample_file)
    ko_merged_df = kegg_reference_df.merge(ko_sample_df, how="right", on="function")
    ko_merged_file_basename = f"{filename_only(ko_sample_file)}_annotated{get_file_extension(ko_sample_file)}"
    ko_merged_file = os.path.join(os.path.dirname(ko_sample_file), ko_merged_file_basename)
    dump_tsv(ko_merged_df, ko_merged_file)

def compose_qiime2_picrust2_results(pipeline_output_dir: str):
    target_dir = os.path.join(os.path.dirname(pipeline_output_dir), "qiime2-picrust2-pipeline-out")
    copying_dicts = deque([
        {
            "source": i,
            "destination": os.path.join(target_dir, "dada2", "qzv", f"dada2_{os.path.basename(i)}"),
        }
        for i in find_by_regex(".+qiime2.+dada2.+\.qzv$", pipeline_output_dir)
    ])
    copying_dicts.extend([
        {
            "source": i,
            "destination": os.path.join(target_dir, "vsearch", "qzv",
                                        f"vsearch_{os.path.basename(i)}"),
        }
        for i in find_by_regex(".+qiime2.+vsearch.+\.qzv$", pipeline_output_dir)
    ])
    for copying_dict in copying_dicts:
        copy(**copying_dict)
    compressing_dicts = deque()
    asv_file = find_by_regex(".+qiime2.+dada2.+ASV_with_taxa.tsv$", pipeline_output_dir)[0]
    compressing_dicts.append({
        "source": asv_file,
        "destination": os.path.join(target_dir, "dada2", "tsv",
                                    f"dada2_{os.path.basename(asv_file)}.7z"),
    })
    otu_normalized_file = \
    find_by_regex(".+qiime2.+vsearch.+OTU_normalized_with_taxa.tsv$", pipeline_output_dir)[0]
    compressing_dicts.append({
        "source": otu_normalized_file,
        "destination": os.path.join(target_dir, "vsearch", "tsv",
                                    f"vsearch_{os.path.basename(otu_normalized_file)}.7z"),
    })
    for picrust2_output_filename in (
            "EC_pred_metagenome_unstrat_described",
            "KO_pred_metagenome_unstrat_described",
            "KO_pred_metagenome_unstrat_described_annotated",
            "pathways_abun_unstrat_described",
    ):
        print(f"Add '{picrust2_output_filename}'")
        picrust2_output_file = \
        find_by_regex(f".+picrust2.+{picrust2_output_filename}\.tsv$", pipeline_output_dir)[0]
        compressing_dicts.append({
            "source": picrust2_output_file,
            "destination": os.path.join(target_dir, "picrust2", "tsv",
                                        f"{os.path.basename(picrust2_output_file)}.7z")
        })
    for compressing_dict in compressing_dicts:
        to_7z(**compressing_dict)


def extract_tables_from_qzv_into_excel(pipeline_output_dir: str):
    qzv_files = find_by_regex(".*\.qzv$", pipeline_output_dir)
    excel_files_dir = os.path.join(pipeline_output_dir, "xlsx")
    for qzv_file in qzv_files:
        qzv_filename = filename_only(qzv_file)
        excel_file = os.path.join(excel_files_dir, f"{qzv_filename}.xlsx")
        extract_tables_from_zip_into_excel(qzv_file, excel_file)


def run(pipeline_output_dir: str, kegg_reference_file: str):
    annotate_qiime2_asv_table(pipeline_output_dir)
    annotate_picrust2_ko_table(pipeline_output_dir, kegg_reference_file)
    compose_qiime2_picrust2_results(pipeline_output_dir)
    extract_tables_from_qzv_into_excel(pipeline_output_dir)


def _parse_args():
    parser = ArgumentParser(
        description="Tool to remove contamination from given FASTA nucleotide sequence file "
                    "based on NCBI contamination report",
        epilog="Note if the FASTA headers contain spaces, only the first piece of split header "
               "will be processed")
    parser.add_argument("-p", "--pipeline", metavar="</path/to/qiime2-picrust2-pipeline>",
                        required=True, help="Input FASTA file")
    parser.add_argument("-k", "--kegg", metavar="</path/to/kegg_denormalized.tsv>", required=True,
                        help="NCBI contamination report file")
    _namespace = parser.parse_args()
    return _namespace.pipeline, _namespace.kegg


if __name__ == '__main__':
    input_pipeline_output_dir, input_kegg_reference_file = _parse_args()
    run(input_pipeline_output_dir, input_kegg_reference_file)
