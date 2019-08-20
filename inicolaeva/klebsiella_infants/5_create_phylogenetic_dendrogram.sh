#!/usr/bin/env bash

echo "Download completed assemblies in GenBank format from NCBI"
DL_DIR=/data1/bio/projects/inicolaeva/klebsiella_infants/ncbi-dl
rm -rf ${DL_DIR}
export IMG=quay.io/biocontainers/ncbi-genome-download:0.2.10--py36_0 && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash -c \
    '
    ping google.com -c 5;
    DL_DIR=/data1/bio/projects/inicolaeva/klebsiella_infants/ncbi-dl;
    mkdir -p $DL_DIR;
    cd $DL_DIR;
    ncbi-genome-download -group bacteria -l complete --genus "Klebsiella pneumoniae" --format genbank --parallel 8 \
    -o $DL_DIR/dl
    chmod -R 777 $DL_DIR
    '
find ${DL_DIR}/dl -name *.gbff.gz -exec gunzip {} \;
mkdir -p ${DL_DIR}/gbff
find ${DL_DIR}/dl -name *gbff -exec mv {} ${DL_DIR}/gbff/ \;


echo "Convert GenBank files to GFF"
rm -rf ${DL_DIR}/gff
export IMG=bioperl/bioperl:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash -c \
    '
    DL_DIR=/data1/bio/projects/inicolaeva/klebsiella_infants/ncbi-dl;
    curl -fsSL https://raw.githubusercontent.com/appris/appris/master/modules/bin/bp_genbank2gff3.pl \
        -o ./bp_genbank2gff3.pl;
    mkdir -p ${DL_DIR}/gff;
    perl bp_genbank2gff3.pl -d ${DL_DIR}/gbff -o ${DL_DIR}/gff;
    chmod -R 777 ${DL_DIR}
    '


echo "Calculate the pan genome"
ROARY_DIR=/data1/bio/projects/inicolaeva/klebsiella_infants/roary
rm -rf ${ROARY_DIR}
mkdir -p ${ROARY_DIR}/gff
find /data1/bio/projects/inicolaeva/klebsiella_infants/pipeline/06_prokka -name *.gff -exec ln -s {} ${ROARY_DIR}/gff/ \;
find /data1/bio/projects/inicolaeva/klebsiella_infants/ncbi-dl/gff -name *.gff -exec ln -s {} ${ROARY_DIR}/gff/ \;
export IMG=sangerpathogens/roary:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash -c \
    '
    ROARY_DIR=/data1/bio/projects/inicolaeva/klebsiella_infants/roary;
    roary -p 32 -f ${ROARY_DIR}/out -e --mafft ${ROARY_DIR}/gff/*.gff;
    chmod -R 777 ${ROARY_DIR}
    '


echo "Get the Newick File"
rm -rf ${ROARY_DIR}/newick
mkdir -p ${ROARY_DIR}/newick
find ${ROARY_DIR}/out -name *.newick | head -n 1 | xargs sed 's|_genomic.gbff||g' - | \
    tee ${ROARY_DIR}/newick/full_raw.newick
chmod -R 777 ${ROARY_DIR}/newick


# Uploaded to iTOOL:
# https://itol.embl.de/tree/53147207280561566215784
# Delete redundant clades, export results to the Newick file 'iTOL_collapsed_tree.newick' and upload it
