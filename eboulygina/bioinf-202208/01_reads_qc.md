# Получение сырых ридов и их контроль качества

## Сырые риды находятся в папке `raw_compressed`

```shell script
ls raw_compressed
```
```text
Saur_S01_L001_R1_001.fastq.gz  Saur_S01_L001_R2_001.fastq.gz
```

## Сырые риды представляют собой сжатые файлы

```shell script
file raw_compressed/Saur_S01_L001_R1_001.fastq.gz
```
```text
raw_compressed/Saur_S01_L001_R1_001.fastq.gz: gzip compressed data, max speed, from NTFS filesystem (NT)
```

```shell script
head -n 2  raw_compressed/Saur_S01_L001_R1_001.fastq.gz
```
```text
<...>
```

##  Распаковка сырых ридов

```shell script
mkdir raw

zcat raw_compressed/Saur_S01_L001_R1_001.fastq.gz > raw/Saur_S01_L001_R1_001.fastq
zcat raw_compressed/Saur_S01_L001_R2_001.fastq.gz > raw/Saur_S01_L001_R2_001.fastq
```

## Распакованные сырые риды представляют собой текстовые файлы, которые имеют определенный формат

```shell script
file raw/Saur_S01_L001_R1_001.fastq
```

```text
raw/Saur_S01_L001_R1_001.fastq: ASCII text
```

```shell script
head -n 5 raw/Saur_S01_L001_R1_001.fastq
```

```text
@M04046:14:000000000-AJC98:1:1101:14884:1380 1:N:0:13
TTTTTTTTTCTTTTTTTTTTTTTTTCTTTTTTTCTCTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTCTTCTTTTCCTTTACTTTTTTTCTTTTCTTTTTTTTTTTTCTTTTTTTCTTTTTTTTTCTTTTTTTTCTTTTTTTTTTCTTTTTTCTTTTTTTTTTTTTTCTTTCTTTTTTCTTTTTTTTCTTTCTTCTTTTTTTTTCTTTCTTTTTTTTTCTTTTTTTTTTTCTCCTTTTTTTTCTTTTTTT
+
1>>>11>>0013BB10AAA//////0DB111//0122B11B>/<01111B/</<>//---<-:-----/0;0;;;C00009000;000-9///////:/;9-----9-/;9/9:9-/;9/:/--9;//:99/9;@/////;-;--9///9;/;////:/9-9;--;--//;////:/99/////9---//9//////9/9/--99///9//99/:---9////;/------/////999/---////;99-
@M04046:14:000000000-AJC98:1:1101:15791:1442 1:N:0:13
```

## Контроль качества сырых ридов осуществляется с помощью `fastqc`

### Развёртывание контейнера `fastqc`

```shell script
# Обновить образ контейнера fastqc версии 0.11.9--hdfd78af_1 с сайта quay.io
docker pull quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1 

docker run \
    --interactive  `# Ждать ввода пользовательской информации в контейнер` \
    --net=host  `# Присоединить контейнер к сети` \
    --rm  `# Не хранить контейнер после выхода` \
    --tty  `# Привязать псевдотерминал к контейнеру` \
    --volume /data1:/data1  `# Смонтировать внешнюю файловую систему хоста в контейнер` \
    --workdir "$(pwd)"  `# Указать, что контейнер должен начать работу в текущей папке` \
    quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1  `# Использовать образ контейнера fastqc версии 0.11.9--hdfd78af_1` \
    bash  `# Перейти в командную оболочку внутри контейнера`
```

### Операции внутри контейнера `fastqc`

```shell script
# Создать конечные папки
mkdir -p \
    fastqc/Saur_S01_L001_R1_001 \
    fastqc/Saur_S01_L001_R2_001

# Запустить программу
fastqc \
    --extract  `# Распаклвать отчеты` \
    --outdir fastqc/Saur_S01_L001_R1_001  `# Указать конечную папку` \
    --threads 20  `# Вести расчеты на 20 потоках` \
    raw/Saur_S01_L001_R1_001.fastq

fastqc \
    --extract \
    --outdir fastqc/Saur_S01_L001_R2_001 \
    --threads 20 \
    raw/Saur_S01_L001_R2_001.fastq

# Найти и далить в пределах папки все файлы с расширением ".zip"
find fastqc/ \
    -name "*.zip" \
    -type f \
    -exec rm "{}" +

# Дать права на чтение/запись остальным пользователям
chmod -fR a+rw fastqc  

# Подсчитать все строки, начинающиеся с `FAIL`
grep -cE '^FAIL' fastqc/Saur_S01_L001_R1_001/Saur_S01_L001_R1_001_fastqc/summary.txt
grep -cE '^FAIL' fastqc/Saur_S01_L001_R2_001/Saur_S01_L001_R2_001_fastqc/summary.txt
```

### Выход из контейнера `fastqc`

```shell script
exit
```
