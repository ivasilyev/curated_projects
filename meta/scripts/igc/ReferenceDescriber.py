#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3

"""

"""
Links:

1. Sequence data for gene catalogs
Integrated non-redundant gene catalog (IGC, nucleotide sequences, fasta)
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.fa.gz

Integrated non-redundant gene catalog (IGC, amino acid sequences, fasta)
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.pep.gz

2. Gene profile, genus profile, KO profile
Gene abundance profile table for 1267 samples (9.9 million genes X 1267 samples)
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/2.ProfileTables/1267sample.gene.relativeAbun.table.gz

Genus profile table for 1267 samples
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/2.ProfileTables/IGC.genus.normalization.ProfileTable.gz

KO profile table for 1267 samples
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/2.ProfileTables/IGC.KO.normalization.ProfileTable.gz

3. Gene annotation summary
IGC annotation and occurrence frequency summary table
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/3.IGC.AnnotationInfo/IGC.annotation_OF.summary.gz

MetaHIT 2010 KEGG annotation
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/6.PreviousGeneCatalog/PGC.gene.KO.list.gz

"""


class ReferenceDescriber:
    name = "IGC"
    description = "Integrated reference catalog of the human gut microbiome"
    documentation = "http://meta.genomics.cn/meta/home"
    # Change the following lines after reference update
    alias = "igc_v2014.03"
    refdata = "/data/reference/IGC/igc_v2014.03/index/igc_v2014.03_refdata.json"
    def export(self):
        print("""
Database alias: {a}
REFDATA linker: {b}
              """.format(a=self.alias, b=self.refdata))
    def parse_refdata(self):
        from meta.scripts.RefDataParser import RefDataParser
        p = RefDataParser(self.refdata)
        return p.get_parsed_list()


class SequenceRetriever:
    def __init__(self):
        import os
        describer = ReferenceDescriber()
        self.alias = describer.alias
        self.reference_dir = "/data/reference/{a}/{b}/".format(a=describer.name, b=self.alias)
        os.makedirs(self.reference_dir, exist_ok=True)
        self.raw_nfasta = "/data/reference/IGC/igc_v2014.03.fasta"
        self.index_dir = "{}index/".format(self.reference_dir)
        self._get_index_guide()
        print("Please replace the following lines in control script:")
        describer.alias = self.alias
        describer.refdata = "{a}{b}.json".format(a=self.index_dir, b=self.alias)
        describer.export()
    def _get_index_guide(self):
        from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
        LaunchGuideLiner.get_index_guide(index_directory=self.index_dir, raw_nfasta_file=self.raw_nfasta)
    @staticmethod
    def _annotate():
        import pandas as pd
        import shutil
        # Retrieved from http://meta.genomics.cn/meta/dataTools
        col_names_list = ["Gene ID",
                          "Gene Name",
                          "Gene Length",
                          "Gene Completeness Status",
                          "Cohort Origin",
                          "Taxonomic Annotation(Phylum Level)",
                          "Taxonomic Annotation(Genus Level)",
                          "KEGG Annotation",
                          "eggNOG Annotation",
                          "Sample Occurence Frequency",
                          "Individual Occurence Frequency",
                          "KEGG Functional Categories",
                          "eggNOG Functional Categories",
                          "Cohort Assembled"]
        annotation_0 = "/data/reference/IGC/igc_v2014.03/index/igc_v2014.03_annotation.tsv"
        shutil.copy2(annotation_0, "{}.bak".format(annotation_0))
        df_0 = pd.read_table(annotation_0, sep='\t', header=0, engine='python')
        df_0["Gene Name"] = df_0.loc[:, "former_id"].apply(lambda x: x.split(" ")[0])
        df_1 = pd.read_table("/data/reference/IGC/IGC.annotation_OF.summary", sep='\t', header='infer', names=col_names_list, engine='python')
        df_3 = pd.merge(df_0, df_1, on="Gene Name", how="outer")
        df_3.to_csv(annotation_0, sep='\t', header=True, index=False)


if __name__ == '__main__':
    retriever = SequenceRetriever()

"""
# Reference indexing (from worker node):

rm -rf /data/reference/IGC/igc_v2014.03/index
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && docker pull $IMG && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG python3 /home/docker/scripts/cook_the_reference.py -i /data/reference/IGC/igc_v2014.03.fasta -n -o /data/reference/IGC/igc_v2014.03/index

Wait until REFDATA file creates

Please replace the following lines in control script:

Database alias: igc_v2014.03
REFDATA linker: /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_refdata.json

"""
