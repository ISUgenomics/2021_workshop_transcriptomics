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

1. Run Salmon

```
....
salmon ...
...
```

## Bee

* reads = `data_bee/reads/*.fastq`
* ref = `data_bee/ref/*.fna.gz`
* gff = `data_bee/ref/*.gff.gz`

```
...
salmon ...
...
```

## Final Gene Counts

Summarize output, how is this different from GSNAP, etc, etc.
