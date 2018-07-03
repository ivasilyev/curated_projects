# Samples preparation
## (Not required) Cut samples to 500 reads
```bash
export S_DIR="/data2/bio/Metagenomes/HG19/Non-mapped_reads"
export D_DIR="/data2/bio/sandbox/201805"
cd ${D_DIR}/..
rm -rf 201805*
mkdir -p ${D_DIR}
for FILE in $(ls -d ${S_DIR}/*.csfasta); do NAME=$(echo ${FILE} | sed "s|${S_DIR}/||g"); head -n 1000 ${FILE} > ${D_DIR}/${NAME}; done
zip -r ${D_DIR}.zip 201805
```

## [On master node] Get reads
```bash
mkdir /data/samples
cd /data/samples
wget -q https://raw.githubusercontent.com/ivasilyev/curated_projects/master/eboulygina/bioinf-201805/201805.zip
unzip 201805.zip
rm -f 201805.zip
```

## Create `sampledata` linker
```bash
cd /data/samples
rm -f 201805.sampledata
export S_DIR="/data/samples/201805"
for FILE in $(ls -d ${S_DIR}/*.csfasta); do printf "$(echo ${FILE} | grep -Po '(?<=201805/).+(?=_no)')\t${FILE}\n" | tee -a 201805.sampledata; done
```

# Reference preparation
## Get reference
```bash
mkdir -p /data/reference/CARD/tmp
cd /data/reference/CARD/tmp
curl -fsSL https://card.mcmaster.ca/latest/data -o card.tar.bz2
tar xvjf card.tar.bz2
```

## Cut reference to 500 sequences
```bash
cd /data/reference/CARD/tmp
head -n 1000 nucleotide_fasta_protein_homolog_model.fasta > ../card_500.fasta
cd ../
rm -rf tmp
```

## [On worker node] Reference indexing
```bash
rm -rf /data/reference/CARD/index
export DOCKER_IMAGE_NAME=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -it ${DOCKER_IMAGE_NAME} python3 /home/docker/scripts/cook_the_reference.py \
-i /data/reference/CARD/card_500.fasta \
-o /data/reference/CARD/index
```

# Job configuration
## Create master config
```bash
mkdir -p /data/charts
cd /data/charts
curl -fsSL https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/config.yaml -o config.yaml
nano config.yaml
```

```text
QUEUE_NAME: card-test-queue
MASTER_CONTAINER_NAME: card-test-master
JOB_NAME: card-test-job
ACTIVE_NODES_NUMBER: 2
THREADS_NUMBER: max
WORKER_CONTAINER_NAME: card-test-worker
SAMPLEDATA: /data/samples/201805.sampledata
REFDATA: /data/reference/CARD/index/card_500.refdata
OUTPUT_MASK: no_hg19_card_500
OUTPUT_DIR: /data/output
```

## Generate charts
```bash
cd /data/charts
curl -fsSL https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/generator.py -o generator.py

export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data --net=host -it ${DOCKER_IMAGE_NAME} python3 /data/charts/generator.py \
-c /data/charts/config.yaml \
-m https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml \
-w https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml \
-o /data/charts
```

# [On master node] Job deployment
## (If not present) Deploy Redis server imperatively
```bash
# View existing pods
kubectl get pods

# Deploy Redis server pod
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml

# Deploy Redis Service
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml
```

## Deploy job imperatively
```bash
cd /data/charts

# Deploy the MASTER chart to create queue
kubectl create -f master.yaml

# Wait until MASTER ("card-test-queue") is completed
kubectl get pods

# If error:
kubectl logs card-test-queue

# Deploy the WORKER chart to create the pipeline job
kubectl create -f worker.yaml
```

# View results
```bash
ls -d /data/output/*

export POD1_NAME=$(kubectl describe pod card-test-job- | grep ^Name: | grep -Po '(?<=card-test-job-).+(?=$)' | head -n 1)
export POD2_NAME=$(kubectl describe pod card-test-job- | grep ^Name: | grep -Po '(?<=card-test-job-).+(?=$)' | tail -n 1)

export POD1_QUEUES_NUM=$(kubectl logs card-test-job-${POD1_NAME} | grep "Loaded full queue on: " | wc -l)
export POD2_QUEUES_NUM=$(kubectl logs card-test-job-${POD2_NAME} | grep "Loaded full queue on: " | wc -l)

printf "\n\n\nTotal files: $(wc -l < /data/samples/201805.sampledata)\nQueue 1 size: ${POD1_QUEUES_NUM}\nQueue 2 size: ${POD2_QUEUES_NUM}\n" && \
printf "Samples with non-zero coverage: $(ls -d /data/output/Statistics/*coverage.txt | wc -l)\n\n\n"
```

# Cleanup
```bash
kubectl delete pod card-test-queue
kubectl delete job card-test-job
```
