# Сборка драфта полногеномной последовательности

```shell script
# Создать конечные папки
mkdir -p spades

# Развернуть контейнер
IMG=quay.io/biocontainers/spades:3.13.1--0 && \
docker pull "${IMG}" && \
docker run --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" sh
```
```shell script
# Найти исполняемый файл
TOOL="$( \
    find /usr/local/ \
        -name "spades.py" \
        -type f \
        2>/dev/null \
    | grep -E 'spades.py$' \
    | head -n 1 \
)"

echo "${TOOL}"

# Запустить программу
python3 "${TOOL}"  `# Запустить через системный интерпретатор Python` \
    --careful  `# Учитывать миссенс- и нонсенс-нуклеотиды` \
    -o spades  `# Выходная папка` \
    --threads "$(grep -c 'processor' /proc/cpuinfo)"  `# Вести расчеты на всех доступных потоках` \
    -1 no_hg38/Saur__unmapped.1.fastq  `# Входные риды` \
    -2 no_hg38/Saur__unmapped.2.fastq \
    2>&1  `# Перенаправить STDERR в STDOUT` \
| tee spades/spades.log  # Добавить логирование

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw spades

# Выйти из контейнера
exit
```

# Контроль качества сборки драфта полногеномной последовательности

```shell script
# Создать конечные папки
mkdir -p quast

# Развернуть контейнер
IMG=quay.io/biocontainers/quast:5.0.2--py36pl5262h30a8e3e_4 && \
docker pull "${IMG}" && \
docker run --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" bash
```

```shell script
# Запустить программу
quast \
    --features quast/sequence.gb  `# Указать референсную последовательность` \
    --gene-finding  `# Включить поиск генов` \
    --min-alignment 200  `# Минимальная длина контига` \
    --no-gzip  `# Не архивировать отчёт` \
    --output-dir quast  `# Выходная папка` \
    --pe1 no_hg38/Saur__unmapped.1.fastq  `# Входные риды для spades` \
    --pe2 no_hg38/Saur__unmapped.2.fastq \
    --plots-format "png"  `# Указать формат изображений` \
    --threads "$(nproc)"  `# Вести расчеты на всех доступных потоках` \
    spades/contigs.fasta  `# Входная сборка` \
|& tee quast/quast.log  # Добавить логирование

# Посмотреть отчёт
cat quast/report.tsv

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw quast

# Выйти из контейнера
exit
```
