# Обрезка ридов с помощью `trimmomatic`

```shell script
# Создать конечные папки
mkdir -p trimmomatic

# Скачать адаптеры
curl \
    "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/tgrigoreva/25_ecoli_genes/adapters.fasta" \
    -o trimmomatic/adapters.fasta

# Развернуть контейнер
IMG="quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4" && \
docker pull "${IMG}" && \
docker run --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" bash
```

```shell script
# Запустить программу
trimmomatic \
    PE  `# Использовать режим для ридов со спаренными концами` \
    -threads 20  `# Вести расчеты на 20 потоках` \
    -phred33  `# Использовать шкалу вероятностей ошибок Phred33` \
    raw/Saur_S01_L001_R1_001.fastq  `# Сырые риды` \
    raw/Saur_S01_L001_R2_001.fastq \
    trimmomatic/Saur_S01_L001_R1_001__trimmed.fastq  `# Обработанные обрезанные риды` \
    trimmomatic/Saur_S01_L001_R1_001__untrimmed.fastq  `# Обработанные необрезанные риды` \
    trimmomatic/Saur_S01_L001_R2_001__trimmed.fastq \
    trimmomatic/Saur_S01_L001_R2_001__untrimmed.fastq \
    ILLUMINACLIP:trimmomatic/adapters.fasta:2:30:10  `# Вырезать указанные адаптеры из ридов` \
    LEADING:3  `# Убрать нуклеотиды с начала, если качество не удовлетворяет` \
    TRAILING:3  `# Убрать нуклеотиды с конца, если качество не удовлетворяет` \
    SLIDINGWINDOW:4:15  `# Убрать нуклеотиды с середины, если качество соседних оснований не удовлетворяет` \
    MINLEN:36  `# Убрать слишком короткие риды` \
| tee "trimmomatic/trimmomatic.log"  # Добавить логирование

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw trimmomatic

# Выйти из контейнера
exit
```

# Удаление адаптеров с помощью `cutadapt`

```shell script
# Создать конечные папки
mkdir -p cutadapt

# Развернуть контейнер
IMG="quay.io/biocontainers/cutadapt:3.5--py37h73a75cf_0" && \
docker pull "${IMG}" && \
docker run -e ADAPTER=AGATCGGAAGAG --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" bash
```

```shell script
# Проверить проброс строки с адаптером
echo "${ADAPTER}"

# Запустить программу
cutadapt \
    --adapter "${ADAPTER}"  `# Удалить прямой адаптер` \
    -A "${ADAPTER}"  `# Удалить компелентарный адаптер` \
    --cores 20  `# Вести расчеты на 20 потоках` \
    --minimum-length 50  `# Убрать слишком короткие риды` \
    --output cutadapt/Saur_S01_L001_R1_001__cut.fastq  `# Выходные файлы` \
    --paired-output cutadapt/Saur_S01_L001_R2_001__cut.fastq \
    trimmomatic/Saur_S01_L001_R1_001__trimmed.fastq `# Входные файлы` \
    trimmomatic/Saur_S01_L001_R2_001__trimmed.fastq \
| tee "cutadapt/cutadapt.log" `# Добавить логирование`

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw cutadapt

# Выйти из контейнера
exit
```
