#!/usr/bin/env bash

# Cleanup
rm -rf /data1/bio/projects/inicolaeva/klebsiella_infants/test/orthomcl
mkdir -p /data1/bio/projects/inicolaeva/klebsiella_infants/test/orthomcl

# Create OrthoMCL sanple data
faa=$(find /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/5_prokka/ -name *.faa -print | sort)
paste <(echo "${faa}") <(echo "${faa}" | grep -Po 'Kleb[0-9]+') | \
    tee /data1/bio/projects/inicolaeva/klebsiella_infants/test/orthomcl/faa.tsv

# Launch the OrthoMCL pipeline container
ping google.com -c 5
export IMG=ivasilyev/orthomcl-mysql:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

# Clean, copy, paste & run the pipeline script
# rm -f /opt/my_tools/pipeline_handler.py && nano /opt/my_tools/pipeline_handler.py
service mysql restart
python3 /opt/my_tools/pipeline_handler.py -i /data1/bio/projects/inicolaeva/klebsiella_infants/test/orthomcl/faa.tsv \
    -o /data1/bio/projects/inicolaeva/klebsiella_infants/test/orthomcl

# One-liner command
export IMG=ivasilyev/orthomcl-mysql:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
bash -c \
    'service mysql restart
     python3 /opt/my_tools/pipeline_handler.py \
        -i /data1/bio/projects/inicolaeva/klebsiella_infants/test/orthomcl/faa.tsv \
        -o /data1/bio/projects/inicolaeva/klebsiella_infants/test/orthomcl'
