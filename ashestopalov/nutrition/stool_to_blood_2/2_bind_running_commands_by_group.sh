
# sudo rm -rf "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline"

mkdir -p -m 0777 "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/scripts/"

cd "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/scripts/"

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/stool_to_blood_2/3_export_variables.sh"

for GROUP_NAME in \
    "Non_obese_adult" \
    "Non_obese_children" \
    "Sport" \
    "NUC_children" \
    "Control" \
    "Obese_adult" \
    "NUC_adult" \
    "Obese_children"
    do
    echo "Analyze '${GROUP_NAME}'"
    bash \
        "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/scripts/3_export_variables.sh" \
        "${GROUP_NAME}" \
        "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/" \
        "/data03/bio/rogachev_human/" \
        "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/sample_data/split/"
    done
