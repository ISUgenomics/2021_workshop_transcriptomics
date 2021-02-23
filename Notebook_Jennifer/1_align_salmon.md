# Salmon - alignment and gene counts

Salmon and Kalisto are both workflow programs. What that means is the indexing, aligning and gene counts are all in one program. (I think...will find out as I run this.)

* [Salmon's official Getting Started Guide](https://combine-lab.github.io/salmon/getting_started/)

## Install Salmon

For now, we'll use `miniconda` to install on Atlas HPC. Salmon is also provided as a precompiled binary (will probably run faster).

1. Create a [salmon_env.yml](bin/salmon_env.yml) miniconda environment file.
 
  ```
  name: salmon_env
  channels:
    - conda-forge
    - bioconda
    - defaults
  dependencies:
    - python=3.8
    - salmon
  ```

2. Create the miniconda environment from Atlas HPC bash

  ```
  conda env create -f salmon_env.yml
  conda env list
  # ... todo add list of envirnments here
  ```
  
3. Make sure `salmon_env` is working

  ```
  conda activate salmon_env
  salmon --help     #<= check usage statement
  ```
  
## Maize

* reads = `data_maize/reads/*_{1,2}.fastq.gz`
* ref = `data_maize/ref/*.fna.gz`
* gff = `data_maize/ref/*.gff.gz`

Oh weird, it doesn't use the full genome, only the transcriptome. I need to fetch the `*.rna.gz` file.

Salmon run is split into two steps: (1) index transcriptome, (2) quantify RNAseq.

1. Index Transcriptome

```
# === Input / Output Variables
REF_RNA=data_maize/ref/*.rna.gz
REF_NAME=b73

salmon index -t ${REF_RNA} -i ${REF_NAME}
```

|flag | value | reason |
|:-:|:-:|:-:|
| -t| data_maize/ref/*.rna.gz| transcriptome|
| -i| b73 | index name, arbitrary genome name |

2. Quantify Reads

```
# === Input / Output Variables
REF_NAME=b73

# === Loop through all read files
for FILE in data_maize/reads/*_1.fastq.gz
do
  READ_NAME=$(basename ${FILE} | sed 's:_1.fastq.gz::g')
  DIR_NAME=$(dirname ${FILE})
  READ_R1=${DIR_NAME}/${READ_NAME}_1.fastq.gz
  READ_R2=${DIR_NAME}/${READ_NAME}_2.fastq.gz
  OUT_COUNTS=${READ_NAME}_genecounts.txt
  
  salmon quant -i ${REF_NAME} \
    -l A \
    -1 ${READ_R1} \
    -2 ${READ_R2} \
    -p 16 \
    --validateMappings \
    -o quants/${OUT_COUNTS}
done
```

|flag | value | reason |
|:-:|:-:|:-:|
| -i| b73 | index name, arbitrary genome name |
| -l| A | |
| -1| *_1.fastq.gz | Left read |
| -2| *_2.fastq.gz | Right read |
| --validateMappings | | a validation step?|
| -o | quants/readname_genecounts.txt | output file|


## Bee

* reads = `data_bee/reads/*.fastq`
* ref = `data_bee/ref/*.fna.gz`
* gff = `data_bee/ref/*.gff.gz`

1. Index Transcriptome

```
# === Input / Output Variables
REF_RNA=data_bee/ref/*.rna.gz
REF_NAME=bombus

salmon index -t ${REF_RNA} -i ${REF_NAME}
```

2. Quantify Reads

```
# === Input / Output Variables
REF_NAME=bombus

# === Loop through all read files
for FILE in data_bee/reads/*.fastq
do
  READ_NAME=$(basename ${FILE} | sed 's:_L002_R1_001.fastq::g')
  DIR_NAME=$(dirname ${FILE})
  READ_R1=${DIR_NAME}/${READ_NAME}_L002_R1_001.fastq
  OUT_COUNTS=${READ_NAME}_genecounts.txt
  
  # Not sure if this will worked with single end, will try...
  salmon quant -i ${REF_NAME} \
    -l A \
    -1 ${READ_R1} \
    -p 16 \
    --validateMappings \
    -o quants/${OUT_COUNTS}
done
```

## Final Gene Counts

Summarize output, how is this different from GSNAP, etc, etc.
