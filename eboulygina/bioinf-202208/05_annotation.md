# Структурная аннотация геномной сборки

```shell script
# Создать конечные папки
mkdir -p prokka

IMG=quay.io/biocontainers/prokka:1.12--pl526_0 && \
docker pull "${IMG}" && \
docker run --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" sh
```

```shell script
prokka \
    --centre UoN  `# Центральный идентификатор` \
    --compliant  `# Режим совместимости с HenBank` \
    --cpu "$(grep -c 'processor' /proc/cpuinfo)"  `# Вести расчеты на всех доступных потоках` \
    --force  `# Перезаписать выходные файлы` \
    --locustag Saur  `# Префикс названия последовательностей оснований` \
    --outdir prokka `# Выходная папка` \
    --prefix Saur  `# Префикс имени файла` \
    --rfam  `# Искать в том числе некодирующие РНК` \
    --genus Staphylococcus  `# Таксономия` \
    --species aureus \
    spades/contigs.fasta  # Входная геномная сборка \
    2>&1  `# Перенаправить STDERR в STDOUT` \
| tee prokka/prokka.log  # Добавить логирование

# Посмотреть статистику аннотации
cat prokka/Saur.txt

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw prokka

# Выйти из контейнера
exit
```

# Функциональная аннотация геномной сборки

## Получение референса

```shell script
# Создать конечные папки
mkdir -p card

# Перейти в конечную папку
cd card

# Скачать и распаковать референсы "на лету"
curl \
    -fsSL \
    "https://card.mcmaster.ca/latest/data" \
| tar jxf  -

# Вернуться в исходную папку
cd ..

# Просмотреть содержимое папки с референсом
ls card
```

## Картирование на референс

```shell script
# Развернуть контейнер
export IMG="quay.io/biocontainers/bowtie2:2.4.4--py37h13ad519_0" && \
docker pull "${IMG}" && \
docker run --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" bash
```

```shell script
# Запустить программу индексации референса
bowtie2-build \
    --threads "$(nproc)"  `# Вести расчеты на всех доступных потоках` \
    card/nucleotide_fasta_protein_homolog_model.fasta  `# Входные референсные последовательности` \
    card/card_mask  `# Маска выходных референсных индексов` \
|& tee card/bowtie2-build.log  # Добавить логирование

# Запустить программу картирования на референс
bowtie2 \
    --gbar 1  `# Убрать результаты с пробелами в районе начала или конца рида` \
    --local  `# Добавить результаты с неполными ридами` \
    --mp 3  `# Убрать результаты со слишком большим числом следующих подряд несовпадений` \
    --threads "$(nproc)"  `# Вести расчеты на всех доступных потоках` \
    -1 no_hg38/Saur__unmapped.1.fastq  `# Входные риды` \
    -2 no_hg38/Saur__unmapped.2.fastq \
    -D 20  `# Число попыток выравнивания каждого участка` \
    -L 3  `# Длина каждого участка выравнивания` \
    -N 1  `# Число попыток выравнивания внутри каждого участка` \
    -R 3  `# Число попыток генерации участков выравнивания` \
    -S card/Saur__mapped.sam  `# Файл с картированными последовательностями ` \
    -x card/card_mask  `# Маска входных референсных индексов` \
|& tee card/bowtie2.log  # Добавить логирование

# Определить долю картировавшихся оснований
grep 'overall alignment rate' card/bowtie2.log

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw card

# Выйти из контейнера
exit
```

## Определение покрытий генов

```shell script
# Развернуть контейнер
export IMG="quay.io/biocontainers/samtools:1.15.1--h1170115_0" && \
docker pull "${IMG}" && \
docker run --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" bash
```

```shell script
# Запустить программу конвертации
samtools view  `# Конвертировать SAM-файл` \
    card/Saur__mapped.sam  `# Входной файл` \
    --bam  `# Формат выходного файла` \
    --uncompressed  `# Не тратить время на сжатие` \
    --threads "$(nproc)"  `# Вести расчеты на всех доступных потоках` \
    2>/dev/null  `# Игнорировать STDERR` \
| samtools sort  `# Сортировать BAM-файл "на лету"` \
    -  `# Использовать STDIN вместо входного файла` \
    -@ "$(nproc)"  `# Вести расчеты на всех доступных потоках` \
    -o card/Saur__sorted.bam  # Выходной файл

# Индексировать BAM-файл
samtools index card/Saur__sorted.bam
```

# Извлечение покрытий

```shell script
# Создать таблицу с заголовком
printf "reference_id\tid_bp\tid_mapped_bp\tid_unmapped_bp\n" \
    > card/samtools_idxstats.tsv

# Добавить к таблице детальную статистику картирования
samtools idxstats \
    card/Saur__sorted.bam \
    >> card/samtools_idxstats.tsv \
    2> card/samtools_idxstats.log

# Извлечь общую статистику картирования
samtools stats \
    card/Saur__sorted.bam \
    > card/samtools_stats.txt \
    2> samtools_stats.log

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw card

# Выйти из контейнера
exit
```

```shell script
# Развернуть контейнер
export IMG="quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6" && \
docker pull "${IMG}" && \
docker run --rm --net=host -it -v /data1:/data1 -w "$(pwd)" "${IMG}" bash
```

```shell script
# Создать таблицу с заголовком
printf "reference_id\tid_coverage_breadth\tid_mapped_bp\tid_bp\tid_mapped_reads_to_id_bp\n" \
    > card/coverages.tsv

# Добавить данные по глубине покрытий к таблице
genomeCoverageBed \
    -ibam `# Использовать BAM в качестве входного формата` \
    card/Saur__sorted.bam \
    2> card/genomeCoverageBed.log \
| grep -Ev '^genome'  `# Убрать строки, начинающиеся с "genome"` \
    >> card/coverages.tsv

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw card

# Выйти из контейнера
exit
```
