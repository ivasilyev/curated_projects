

export IMG="ivasilyev/curated_projects:latest"
docker pull "${IMG}"
docker run \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./meta/scripts/process_qiime2_picrust2_results.py \
                --pipeline "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline" \
                --kegg "/data/reference/KEGG/kegg_v2024-05-11/kegg_v2024-05-11_denormalized.tsv"
        '
