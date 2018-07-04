#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

import subprocess
import os
import requests
import bs4
import lxml
import re
import multiprocessing
import pandas as pd
import yaml

# E-utilities Quick Start: https://www.ncbi.nlm.nih.gov/books/NBK25500/
# Valid values of &retmode and &rettype for EFetch: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly


class Utilities:
    @staticmethod
    def get_page(url):
        header = "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US) AppleWebKit/525.19 (KHTML, like Gecko) Chrome/1.0.154.53 Safari/525.19"
        return requests.get(url, headers={'User-Agent': header}).content
    @staticmethod
    def remove_empty_values(input_list):
        output_list = []
        if input_list is not None:
            for i in input_list:
                if i is not None:
                    try:
                        if len(i) > 0:
                            output_list.append(i)
                    except TypeError:
                        continue
        return output_list
    @staticmethod
    def dict2pd_series(dictionary):
        output = pd.Series()
        for key in dictionary:
            output.at[key] = dictionary[key]
        return output
    @staticmethod
    def multi_core_queue(function_to_parallelize, queue):
        pool = multiprocessing.Pool()
        output = pool.map(function_to_parallelize, queue)
        pool.close()
        pool.join()
        return output


class FASTA:
    """
    This class is an attempt to apply NCBI standards to single FASTA.
    Consumes one header followed by sequence.
    """
    def __init__(self, single_fasta):
        self._body = re.sub("\n+", "\n", single_fasta.replace('\r', ''))
        try:
            self.header = re.findall("^>(.+)", self._body)[0].strip()
            self.sequence = "\n".join(self.chunk_string(re.sub("[^A-Za-z]", "", self._body.replace(self.header, "")), 70)).strip().upper()
        except IndexError:
            raise ValueError("Cannot parse the header for sequence: {}".format(self._body))
        # Nucleotide sequence has only AT(U)GC letters. However, it may be also protein FASTA.
    @staticmethod
    def chunk_string(string, length):
        return [string[0 + i:length + i] for i in range(0, len(string), length)]
    def __len__(self):
        return len(self.sequence)
    def to_dict(self):
        return {self.header: self.sequence}
    def to_str(self):
        return "\n".join([">" + self.header, self.sequence])
    def set_header(self, header):
        if header:
            if len(header) > 0:
                self.header = header


class FASTAParser:
    """
    This class parses FASTA file read as single string
    """
    def __init__(self, fastas_string):
        self._fastas_string = fastas_string
        self._raw_fastas_list = [">{}".format(j) if not j.startswith(">") else j for j in [i.strip() for i in re.split("\n>", self._fastas_string)]]
        self._parsed_fastas_list = Utilities.remove_empty_values([FASTA(i) for i in self._raw_fastas_list])
    def get_full_length(self):
        return sum([len(i) for i in self._parsed_fastas_list])
    def get_full_sequence(self):
        return "\n".join([i.to_str() for i in self._parsed_fastas_list])
    def get_fastas_list(self):
        return self._parsed_fastas_list
    def get_headers_list(self):
        return [i.header for i in self._parsed_fastas_list]


class TadbFastaHeaderAnnotator:
    """
    This class parses TADB feature page (http://202.120.12.135/TADB2/feature_page.php?TAs_id=<feature_id>)
    Consumes FASTA header
    Returns Pandas Series supplied with parsed details of nucleotide and protein FASTA headers
    """
    def __init__(self, header):
        self._header = header
        self._annotation_dict = self._parse_nfasta_header(self._header)
        self._url = "http://202.120.12.135/TADB2/feature_page.php?TAs_id={}".format(self._annotation_dict["tadb_id"])
        self._soup = bs4.BeautifulSoup(Utilities.get_page(self._url), "lxml")
        # Parse 'TA locus' tab
        self._ta_locus_dict = {k: "" for k in ['TA ID', 'Type', 'Replicon', 'Strain', 'Location', 'Classification (toxin-antitoxin)', 'Toxin', 'Antitoxin', "Regulator"]}
        for tr_soup in self._soup.find_all("div", class_="panes")[0].find_all("table")[0].find_all("tr"):
            try:
                self._ta_locus_dict[tr_soup.find_all("td")[0].text.strip()] = tr_soup.find_all("td")[1].text.strip()
            except IndexError:
                continue
        # Parse 'Sequence' tab
        sequence_url = "http://202.120.12.135/TADB2/feature_page.php?TAs_id={}#Sequence".format(self._annotation_dict["tadb_id"])
        web_soup = self._soup.find_all(face="Courier", size="-1")
        web_fasta = ""
        try:
            if self._annotation_dict["category"] == "Toxin":
                web_fasta = web_soup[0].get_text("\n")
            elif self._annotation_dict["category"] == "Antitoxin":
                web_fasta = web_soup[2].get_text("\n")
            elif self._annotation_dict["category"] == "Regulator":
                try:
                    web_fasta = web_soup[4].get_text("\n")
                except IndexError:
                    print("The page {} is categorized as '{}', but does not supply the required sequence. Falling back to 'Antitoxin' category.".format(sequence_url, self._annotation_dict["category"]))
                    web_fasta = web_soup[2].get_text("\n")
        except IndexError:
            print(web_soup)
            print(len(web_soup))
            raise ValueError("Cannot parse the sequences. Please check the page: {}".format(sequence_url))
        self._protein_fasta = FASTA(web_fasta)
        self._protein_header = "{} {}".format(re.split(" ", self._header)[0], re.sub("^(.* )gi|^(.* )ref", "", self._protein_fasta.header))
        self._protein_header_dict = self._parse_pfasta_header(self._protein_header)
        self._annotation_dict.update({k: self._protein_header_dict[k] for k in self._protein_header_dict if k in ["genpept_id", "product"]})
        self._annotation_dict.update({"type": self._ta_locus_dict["Type"],
                                      "t/at_classification": self._ta_locus_dict["Classification (toxin-antitoxin)"],
                                      "toxin_properties": self._ta_locus_dict["Toxin"],
                                      "antitoxin_properties": self._ta_locus_dict["Antitoxin"],
                                      "regulator_properties": self._ta_locus_dict["Regulator"]})
        self._annotation_dict.update({"raw_nfasta_header": self._header, "raw_pfasta_header": self._protein_header})
    @staticmethod
    def _parse_nfasta_header(header):
        for tag in re.findall("\[(.+)\]", header):
            header = header.replace("[{}]".format(tag), "[{}]".format(tag.strip()))
        tadb_chunk = re.split(" ", header)[0].upper().strip()
        ids_chunk = re.split(" ", header)[1]
        strain_chunk = " ".join(re.split(" ", header)[2:-1])
        gene_name_chunk = re.split(" ", header)[-1]
        output_dict = {"tadb_id": re.findall("^TADB\|[A-Z]{1,2}([0-9]+)", tadb_chunk)[0],
                       "category": re.findall("^TADB\|([A-Z]{1,2})[0-9]+", tadb_chunk)[0]}
        if output_dict["category"] == "T":
            output_dict["category"] = "Toxin"
        elif output_dict["category"] == "AT":
            output_dict["category"] = "Antitoxin"
        elif output_dict["category"] == "RE":
            output_dict["category"] = "Regulator"
        else:
            raise ValueError("Cannot define the header's category: {}".format(header))
        output_dict.update({"geninfo_id": re.findall("gi\|([0-9]+)\|", ids_chunk),
                            "genbank_id": re.findall("ref\|(.+)\|", ids_chunk),
                            "locus": re.findall("\|:[A-Za-z]*([0-9]+\-[0-9]+)", ids_chunk),
                            "is_antisense_strand": re.findall("\|:([A-Za-z]*)[0-9]+\-[0-9]+", ids_chunk),
                            "strain": strain_chunk,
                            "gene_name": re.findall("\[(.+)\]$", gene_name_chunk)})
        output_dict = {k: "".join(output_dict[k]) for k in output_dict}
        output_dict["is_antisense_strand"] = output_dict["is_antisense_strand"].lower() == "c"
        output_dict["strain"] = re.sub(" +", " ", output_dict["strain"].replace(output_dict["genbank_id"], "").strip())
        return output_dict
    @staticmethod
    def _parse_pfasta_header(header):
        d = TadbFastaHeaderAnnotator._parse_nfasta_header(header)
        replacing_dict = {"genbank_id": "genpept_id", "strain": "product", "gene_name": "protein_name"}
        for k in replacing_dict:
            d[replacing_dict[k]] = d.pop(k)
        return d
    def get_annotation(self):
        if self._annotation_dict:
            return Utilities.dict2pd_series(self._annotation_dict)
        else:
            return


class TadbFastaAnnotator:
    """
    This class parses TADB nucleotide FASTA
    Consumes original FASTA file URL (http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/<feature_id>)
    Returns FASTA headers
    """
    def __init__(self, url):
        self._url = url
        # self_url = 'http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_AT_exp.fas'
        self._whole_parser = FASTAParser(Utilities.get_page(self._url).decode('utf-8'))
        self._headers_queue = Utilities.remove_empty_values(self._whole_parser.get_headers_list())
        self._mp_annotations_list = Utilities.remove_empty_values(Utilities.multi_core_queue(self._mp_annotate_header, self._headers_queue))
    @staticmethod
    def _mp_annotate_header(header):
        try:
            annotator = TadbFastaHeaderAnnotator(header)
            return annotator.get_annotation()
        except Exception as err:
            print("Failed to process header: {}".format(header))
            raise err
    def get_annotation_table(self):
        return pd.concat(self._mp_annotations_list, axis=1).transpose().set_index("tadb_id")
    def get_fastas_list(self):
        return self._whole_parser.get_fastas_list()
    """
    ann = TadbFastaAnnotator("http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_AT_exp.fas")
    ann.get_annotation_table()

    Collected exceptions: 
    
    The page http://202.120.12.135/TADB2/feature_page.php?TAs_id=6104#Sequence is categorized as 'Regulator', but does not supply the required sequence. Falling back to 'Antitoxin' category.
    Failed to process header: TADB|T6272 gi|15600771|ref|NC_002506.1|:333953-334219 NC_002506.1 Vibrio cholerae O1 biovar El Tor str. N16961 chromosome II, complete sequence [VCA0332]
    Failed to process header: TADB|T230 gi|16124256|ref|NC_002696.2|:c2952832-2951498 NC_002696.2 Caulobacter crescentus CB15 chromosome, complete genome [CC_2735]
    Failed to process header: TADB|Tgi|81230333|ref|NC_007595.1|:21934-22212 NC_007595.1 Synechococcus elongatus PCC 7942 plasmid 1, complete sequence [Synpcc7942_B2631]
    Failed to process header: TADB|Tgi|94987631|ref|NC_008021.1|:453674-453979 NC_008021.1 Streptococcus pyogenes MGAS9429, complete genome [MGAS9429_Spy0456]
    Failed to process header: TADB|Tgi|121999251|ref|NC_008790.1|:17771-17974 NC_008790.1 Campylobacter jejuni subsp. jejuni 81-176 plasmid pTet, partial sequence [CJJ81176_pTet0022]
    Failed to process header: TADB|Tgi|134047083|ref|NC_009140.1|:c15404-15054 NC_009140.1 Salmonella enterica subsp. enterica serovar Newport str. SL254 plasmid pSN254, complete sequence [SNSL254_pSN254_0025]
    Failed to process header: TADB|Tgi|161789063|ref|NC_010115.1|:c15113-14748 NC_010115.1 Coxiella burnetii RSA 331 plasmid QpH1, complete sequence [COXBURSA331_0025]
    Failed to process header: TADB|Tgi|170076470|ref|NC_010474.1|:c40774-40448 NC_010474.1 Synechococcus sp. PCC 7002 plasmid pAQ7, complete sequence [SYNPCC7002_G0034]
    Failed to process header: TADB|Tgi|387619774|ref|NC_017638.1|:2743827-2744900 NC_017638.1 Escherichia coli DH1 (ME8569) DNA, complete genome [ECDH1ME8569_2549]
    Failed to process header: TADB|Tgi|190895317|ref|NC_010998.1|:43895-44179 NC_010998.1 Rhizobium etli CIAT 652 plasmid pA, complete sequence [RHECIAT_PA0000034]
    Failed to process header: TADB|T6179 gi|78221228|ref|NC_007517.1|:354256-354594 NC_007517.1 Geobacter metallireducens GS-15, complete genome [Gmet_0316]
    Failed to process header: TADB|Tgi|190889639|ref|NC_010994.1|:c296889-296590 NC_010994.1 Rhizobium etli CIAT 652, complete genome [RHECIAT_CH0000281]
    Failed to process header: TADB|T6179 gi|78221228|ref|NC_007517.1|:354256-354594 NC_007517.1 Geobacter metallireducens GS-15, complete genome [Gmet_0316]
    Failed to annotate the remote FASTA link: http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_T.fas
    ...
    """


class TADBDownloadPageParser:
    """
    This class parses TADB download page (http://202.120.12.135/TADB2/download.html) and performs advanced annotation
    """
    def __init__(self):
        self._url = "http://202.120.12.135/TADB2/download.html"
        self._soup = bs4.BeautifulSoup(Utilities.get_page(self._url), "lxml")
        self._table = self._soup.find_all(id="bodyPanel")[0]
        # The experimentally validated TA loci have 'bgcolor="#EEFFEE"'
        # All the in silico predicted and experimentally validated TA loci have 'bgcolor="#d0e0e3"'
        self._nucleotide_fasta_links_list = ["http://202.120.12.135/TADB2/{}".format(i["href"]) for i in self._table.find_all(href=re.compile("nucleotide"))]
        annotators_list = []
        for link in self._nucleotide_fasta_links_list:
            try:
                annotators_list.append(TadbFastaAnnotator(link))
            except Exception as err:
                print("Failed to annotate the remote FASTA link: {}".format(link))
                raise err


class NewTadbDownloadPageParser:
    """
    This class parses TADB download page (http://202.120.12.135/TADB2/download.html) and performs simplified annotation
    """
    def __init__(self):
        self._url = "http://202.120.12.135/TADB2/download.html"
        self._soup = bs4.BeautifulSoup(Utilities.get_page(self._url), "lxml")
        self._table_rows_soups_list = Utilities.remove_empty_values([i.find_all("a") for i in self._soup.find_all(id="bodyPanel")[0].find_all("tr")])
        self._table_rows_dicts_list = []
        # Prepare dicts '{<FASTA category>: <FASTA URL>}'
        for table_row_soup in self._table_rows_soups_list:
            self._table_rows_dicts_list.append({re.sub("[^A-Za-z ]", "", i.text).strip(): "http://202.120.12.135/TADB2/{}".format(i["href"]) for i in table_row_soup})
        # Parse FASTA URLs to FASTAParser objects
        self._fastas_dicts_list = [{k: FASTAParser(Utilities.get_page(i[k]).decode('utf-8')) for k in i} for i in self._table_rows_dicts_list]
        # Combine whole protein FASTA headers list
        pfastas_raw_headers_combined_list = []
        for fastas_dict in self._fastas_dicts_list:
            pfastas_parser = fastas_dict.get("Amino acid")
            if pfastas_parser:
                pfastas_raw_headers_combined_list.extend(pfastas_parser.get_headers_list())
        # Try to perform rough matches between nFASTA and pFASTA headers
        matching_dicts_list = []
        for fastas_dict in self._fastas_dicts_list:
            nfastas_parser = fastas_dict.get("Nucleotide")
            if nfastas_parser:
                nfastas_list = nfastas_parser.get_fastas_list()
                _mp_search_pfasta_by_nfasta_queue = [{"nfasta": i, "pfasta_headers_list": pfastas_raw_headers_combined_list} for i in nfastas_list]
                matching_dicts_list.extend(Utilities.multi_core_queue(self._mp_search_pfasta_by_nfasta, _mp_search_pfasta_by_nfasta_queue))
        # Process matched and unmatched events
        self._annotation_table = pd.concat([i["series"] for i in matching_dicts_list], axis=1).transpose().set_index("nfasta_header")
        self._fastas_list = [i["nFASTA"] for i in matching_dicts_list]
    @staticmethod
    def _mp_search_pfasta_by_nfasta(input_dict):
        nfasta_header = input_dict["nfasta"].header
        matching_chunk = re.split(" ", nfasta_header)[0]
        matching_pfasta_header = ""
        gi_id = re.findall("gi\|([0-9]+)\|", nfasta_header)
        nfasta_ref_id = re.findall("ref\|([0-9]+)\|", nfasta_header)
        for pfasta_header in input_dict["pfasta_headers_list"]:
            if pfasta_header.startswith(matching_chunk):
                matching_pfasta_header = pfasta_header
                break
            elif len(gi_id) > 0:
                if "gi|{}|".format(gi_id[0]) in pfasta_header:
                    matching_pfasta_header = pfasta_header
            elif len(nfasta_ref_id) > 0:
                pfasta_ref_id = re.findall("ref\|.([0-9]+)\|", pfasta_header)
                if len(pfasta_ref_id) > 0:
                    if nfasta_ref_id[0] == pfasta_ref_id[0]:
                        matching_pfasta_header = pfasta_header
        return {"series": Utilities.dict2pd_series({"nfasta_header": nfasta_header, "pfasta_header": matching_pfasta_header}),
                "nFASTA": input_dict["nfasta"]}
    def get_annotation(self):
        return self._annotation_table
    def get_sequence(self):
        return "\n".join([i.to_str() for i in self._fastas_list]) + "\n"
    def dump_annotation(self, file):
        df = self.get_annotation()
        df.to_csv(file, sep='\t', index=True, header=True)
    def dump_sequence(self, file):
        with open(file, 'w') as f:
            f.write(self.get_sequence())
    def update_annotation(self, file):
        annotation_subdf = self.get_annotation()
        annotation_subdf.index.names = ["former_id"]
        annotation_df = pd.read_table(file, sep='\t', header=0, engine='python')
        annotation_df.to_csv("{}.bak".format(file), sep='\t', index=False, header=True)
        out_df = pd.concat([annotation_df.set_index("former_id"), annotation_subdf], axis=1)
        out_df = out_df.reset_index().set_index("reference_id").sort_index(ascending=True)
        out_df.to_csv(file, sep='\t', index=True, header=True)


if __name__ == '__main__':
    parser = NewTadbDownloadPageParser()

    referenceDir = "/data/reference/TADB/"
    referenceDBName = "tadb_v2.0"
    referenceSequenceFile = "{referenceDir}{referenceDBName}.fasta".format(referenceDir=referenceDir, referenceDBName=referenceDBName)
    os.makedirs(referenceDir, exist_ok=True)

    parser.dump_sequence(referenceSequenceFile)

    print("""
    # Reference indexing
    rm -rf {referenceDir}index
    docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest \
    python3 /home/docker/scripts/cook_the_reference.py \
    -i {referenceSequenceFile} \
    -o {referenceDir}index
    """.format(referenceDir=referenceDir, referenceSequenceFile=referenceSequenceFile))

    """
    # Reference indexing
    rm -rf /data/reference/TADB/index
    docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest python3 /home/docker/scripts/cook_the_reference.py -i /data/reference/TADB/tadb_v2.0.fasta -o /data/reference/TADB/index
    """

    parser.update_annotation("/data/reference/TADB/index/tadb_v2.0_annotation.txt")

    ownerName = "dsafina"
    projectName = "hp_checkpoints/TADB"
    outputDir = "/data1/bio/projects/{a}/{b}/".format(a=ownerName, b=projectName)
    chartsDir = "{}charts/".format(outputDir)
    deployName = "dsafina-bwt-tadb"

    cfgDict = {"QUEUE_NAME": "{}-queue".format(deployName),
               "MASTER_CONTAINER_NAME": "{}-master".format(deployName),
               "JOB_NAME": "{}-job".format(deployName),
               "ACTIVE_NODES_NUMBER": 9,
               "THREADS_NUMBER": "half",
               "WORKER_CONTAINER_NAME": "{}-worker".format(deployName),
               "SAMPLEDATA": "/data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata",
               "REFDATA": "{referenceDir}index/{referenceDBName}.refdata".format(referenceDir=referenceDir, referenceDBName=referenceDBName),
               "OUTPUT_MASK": "no_hg19",
               "OUTPUT_DIR": "/data2/bio/Metagenomes/Toxins/TADB"}

    cfgFileName = "{}config.yaml".format(chartsDir)

    subprocess.getoutput("mkdir -p {}".format(chartsDir))
    with open(cfgFileName, 'w') as cfgFile:
        yaml.dump(cfgDict, cfgFile, default_flow_style=False, explicit_start=False)

    genFileName = chartsDir + "generator.py"
    subprocess.getoutput("curl -fsSL https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/generator.py -o {}".format(genFileName))
    subprocess.getoutput("python3 {genFileName} -c {cfgFileName} -m https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml -w https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml -o {chartsDir}".format(genFileName=genFileName, cfgFileName=cfgFileName, chartsDir=chartsDir))
    os.remove(genFileName)

    print("""
    # Copy the directory '{chartsDir}' into '{ownerName}/{projectName}' directory and push updates
    
    # Look for Redis pod & service:
    kubectl get pods --show-all
    
    # Deploy if not present:
    kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && \
    kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml
    
    # Pipeline launch
    # Deploy the MASTER chart to create queue
    kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/{ownerName}/{projectName}/charts/master.yaml
    
    # Wait until master finish and deploy the WORKER chart to create the pipeline job
    kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/{ownerName}/{projectName}/charts/worker.yaml
    
    # View active nodes
    kubectl describe pod {JOB_NAME}- | grep Node:
    
    # View progress (from WORKER node)
    echo && echo PROCESSED $(ls -d {OUTPUT_DIR}/Statistics/*coverage.txt | wc -l) OF $(cat {sampleDataFileName} | wc -l)
    
    # Look for some pod (from MASTER node)
    kubectl describe pod <NAME>
    
    # Cleanup
    kubectl delete pod {QUEUE_NAME}
    kubectl delete job {JOB_NAME}
    
    # Checkout (from WORKER node)
    docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 \
    /home/docker/scripts/verify_coverages.py -s {sampleDataFileName} \
    -r {REFDATA} \
    -m {OUTPUT_MASK}_{referenceDBName} -d -o {OUTPUT_DIR}
    
    """.format(ownerName=ownerName,
               chartsDir=chartsDir,
               projectName=projectName,
               JOB_NAME=cfgDict["JOB_NAME"],
               OUTPUT_DIR=cfgDict["OUTPUT_DIR"],
               QUEUE_NAME=cfgDict["QUEUE_NAME"],
               sampleDataFileName=cfgDict["SAMPLEDATA"],
               REFDATA=cfgDict["REFDATA"],
               OUTPUT_MASK=cfgDict["OUTPUT_MASK"],
               referenceDBName=referenceDBName))

    """
    # Copy the directory '/data1/bio/projects/dsafina/hp_checkpoints/TADB/charts/' into 'dsafina/hp_checkpoints/TADB' directory and push updates

    # Look for Redis pod & service:
    kubectl get pods --show-all
    
    # Deploy if not present:
    kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml
    
    # Pipeline launch
    # Deploy the MASTER chart to create queue
    kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/dsafina/hp_checkpoints/TADB/charts/master.yaml
    
    # Wait until master finish and deploy the WORKER chart to create the pipeline job
    kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/dsafina/hp_checkpoints/TADB/charts/worker.yaml
    
    # View active nodes
    kubectl describe pod dsafina-bwt-tadb-job- | grep Node:
    
    # View progress (from WORKER node)
    echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/Toxins/TADB/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata | wc -l)
    
    # Look for some pod (from MASTER node)
    kubectl describe pod <NAME>
    
    # Cleanup
    kubectl delete pod dsafina-bwt-tadb-queue
    kubectl delete job dsafina-bwt-tadb-job
    
    # Checkout (from WORKER node)
    docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 /home/docker/scripts/verify_coverages.py -s /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata -r /data/reference/TADB/index/tadb_v2.0.refdata -m no_hg19_tadb_v2.0 -d -o /data2/bio/Metagenomes/Toxins/TADB

    """


class GroupDataPreparer:
    def __init__(self, raw_groupdata_file, processed_prefix, processed_suffix):
        self._df = pd.read_table(filepath_or_buffer=raw_groupdata_file,
                                 sep='\t',
                                 header='infer',
                                 names=["sample_name", "group_name"],
                                 engine='python')
        self._df["file_name"] = self._df.loc[:, "sample_name"].map(lambda x: processed_prefix + x + processed_suffix)
        self._df["file_exists"] = self._df.loc[:, "file_name"].map(lambda x: os.path.isfile(x))
        self._df = self._df.loc[self._df.loc[:, "file_exists"] == True]
        if sum(self._df.loc[:, "file_exists"].values.tolist()) == 0:
            raise ValueError("Cannot find files by the mask: '{a}<sample name>{b}'".format(a=processed_prefix, b=processed_suffix))
    def get_groupdata(self):
        return self._df.loc[:, ["file_name", "group_name"]]
    def export_groupdata(self, output_file):
        self.get_groupdata().to_csv(path_or_buf=output_file, sep='\t', header=False, index=False)


class RemoteScript:
    def __init__(self, url, output_dir):
        self._url = url
        if not output_dir:
            output_dir = os.getcwd()
        output_dir = (output_dir + "/", output_dir)[output_dir.endswith("/")]
        script_name = url.split("/")[-1]
        self.file = output_dir + script_name
        subprocess.getoutput("curl -fsSL {a} -o {b}".format(a=self._url, b=self.file))


class PivotTableAnnotator:
    def __init__(self, pivot_file, annotation_file):
        pivot_df = pd.read_table(pivot_file, sep='\t', header=0, engine='python')
        annotation_df = pd.read_table(annotation_file, sep='\t', header=0, engine='python')
        self.annotated_df = pd.merge(annotation_df, pivot_df, on="reference_id", how='outer')
    def export_annotated_df(self, output_file):
        self.annotated_df.to_csv(path_or_buf=output_file, sep='\t', header=True, index=False)



def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


groupdataFile = "/data1/bio/projects/dsafina/HP-checkpoints_SRR.groupdata"
ownerName = "dsafina"
projectName = "hp_checkpoints"
referenceDBName = "TADB"
projectsDir = "/data1/bio/projects/"
outputDir = "{pdir}{o}/{pname}/{r}/".format(pdir=projectsDir, o=ownerName, pname=projectName, r=referenceDBName)
os.makedirs(path=outputDir, exist_ok=True)
preparer = GroupDataPreparer(raw_groupdata_file=groupdataFile,
                             processed_prefix="/data2/bio/Metagenomes/Toxins/TADB/Statistics/",
                             processed_suffix="_no_hg19_tadb_v2.0_coverage.tsv")
processed_group_data = "{a}{b}.groupdata".format(a=outputDir, b=filename_only(groupdataFile))
preparer.export_groupdata(processed_group_data)

g2s_script = RemoteScript("https://raw.githubusercontent.com/ivasilyev/statistical_tools/master/groupdata2statistics.py", None)
subprocess.getoutput("python3 {a} -g {g} -i {i} -v {v} -o {o}".format(a=g2s_script.file,
                                                                      g=processed_group_data,
                                                                      i="reference_id",
                                                                      v="id_total_relative_abundance",
                                                                      o=outputDir + filename_only(g2s_script.file)))

annotator = PivotTableAnnotator(pivot_file="/data1/bio/projects/dsafina/hp_checkpoints/TADB/groupdata2statistics/C_I_II_III_IV_srr_total_dataframe.tsv", annotation_file="/data/reference/TADB/index/tadb_v2.0_annotation.tsv")
annotator.export_annotated_df("/data1/bio/projects/dsafina/hp_checkpoints/TADB/groupdata2statistics/C_I_II_III_IV_srr_total_dataframe_annotated.tsv")
