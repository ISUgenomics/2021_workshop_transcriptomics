# Methods

Written summary of methods performed in this repo. A lot of the steps described here were copied from Jennifer's notebook - thanks Jennifer!! I added more details to help myself in the future since I'm a newbie ^_^

## Raw data
* **Maize data:** https://www.ebi.ac.uk/ena/browser/view/PRJNA260793
* **Maize reference (*Zea mays* B73):** ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.*gz
* **Bee data**
* **Bee reference (*Bombia impatiens*):** https://hymenoptera.elsiklab.missouri.edu/genome_fasta<br /> Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa

## Data transfer to HPC (Atlas dtn node)
* `ssh username@Atlas-dtn.hpc.msstate.edu`
https://www.hpc.msstate.edu/computing/atlas/
  * dtn node: faster speed at transferring data onto the HPC

### Maize
`/projectdirectory/mydirectory/rnaseq/maize`
1. Fetch maize reference genome:
```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.*gz
```

2. Slurm script for fetching maize sequences to Atlas. See Jennifer's `atlas_maizedata.slurm` for template
<details><summary>See Jennifer's `atlas_maizedata.slurm` for template</summary>

  ```
#!/bin/bash
#SBATCH --job-name=Maize                             # name of the job submitted
#SBATCH -p service                                   # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 2                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 24:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru
#Enter commands here:
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573506/SRR1573506_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573506/SRR1573506_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573507/SRR1573507_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573507/SRR1573507_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573508/SRR1573508_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573508/SRR1573508_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573509/SRR1573509_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573509/SRR1573509_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573510/SRR1573510_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573510/SRR1573510_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573511/SRR1573511_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573511/SRR1573511_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573512/SRR1573512_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573512/SRR1573512_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573513/SRR1573513_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573513/SRR1573513_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573514/SRR1573514_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573514/SRR1573514_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573515/SRR1573515_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573515/SRR1573515_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573516/SRR1573516_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573516/SRR1573516_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573517/SRR1573517_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573517/SRR1573517_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573518/SRR1573518_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573518/SRR1573518_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573519/SRR1573519_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573519/SRR1573519_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573520/SRR1573520_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573520/SRR1573520_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573521/SRR1573521_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573521/SRR1573521_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573522/SRR1573522_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573522/SRR1573522_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573523/SRR1573523_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573523/SRR1573523_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573524/SRR1573524_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573524/SRR1573524_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573525/SRR1573525_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573525/SRR1573525_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573526/SRR1573526_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573526/SRR1573526_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573527/SRR1573527_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573527/SRR1573527_2.fastq.gz
#End of file
```
  `Submitted batch job 127697`
  * Fetching sequences with `wget`: https://github.com/ISUgenomics/2021_workshop_transcriptomics/blob/main/00_Files.md
  </details>

### Bee
`/projectdirectory/mydirectory/rnaseq/bee`

1. Download all data locally in `bee/`

2. Copy fastq files from raw data directories to `bee/`:
```
find . -name *.fastq -exec cp '{}' "./raw/" ";"`
```

3. Copy .sam and Bombus_impatiens_* files from mapping directories to `mapping/`:
```
find . -name *.sam -exec cp '{}' "./mapping/" ";"`
find . -name Bombus_impatiens_* -exec cp '{}' "./mapping/" ";"
```
  * Bee reference genome in here: `Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa`

4. Tar raw and mapping directories together
```
tar -zcvf bee.tar.gz bee/
```

5. Uploaded `bee.tar.gz` to Atlas. Takes about 1.5h for 22GB.

## File and software tool locations on Atlas
```
ProjectDirectory/
  |_MyDirectory/
      |_bee/
          |_bee.tar.gz
          |_logs/
          |_mapping/
          |_meta/
          |_raw_data/
          |_reference_genome_bee/
          |_results/
              |_gsnap/
          |_scripts/
      |_maize/
          |_gsnap/
              |_b73_reference_gsnap/
          |_logs/
          |_metadata/
          |_raw_data/
          |_reference_genome/
          |_results/    
          |_scripts/
      |_dot_files
          |_miniconda3/
          |_Miniconda3-latest-Linux-x86_64.sh
          |_software/
              |_bbmap/
              |_BBMap_38.86.tar
              |_gmap # linked to executable>
              |_gmap-2020-12-17/
                  |_bin/ # executables in here
                      |_gmap_build
                      |_gmap
              |_gmapdb/
                  |_b73/
              |_gmap-gsnap-2020-12-17.tar
```
```
home/
    |_rnaseq # linked to projectdirectory/mydirectory/
    |_inbox/ #linked to projectdirectory/mydirectory/
    |_outbox/ #linked to projectdirectory/mydirectory/
    |_miniconda3/ #linked to projectdirectory/mydirectory/dot_files/
    |_software # linked to projectdirectory/mydirectory/dot_files/
        |_bbmap/
        |_BBMap_38.86.tar
        |_gmap # symbolic link to gmap-2020-12-17/src/gmap.avx2>
        |_gmap-2020-12-17/
        |_gmap-gsnap-2020-12-17.tar
```
* dotfiles (`.singularity, .conda`) are usually invisible folders that get large as you install conda packages, or singularity images. These can eat up your home folder ~5GB memory limit if they're not softlinked

## Install gsnap on Atlas -- giving up because too difficult to locally install other tools that gsnap needs

<details><summary>Notes from a local **GSNAP** install</summary>

1. Fetch gsnap software to Atlas in project directory.
```
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-12-17.tar.gz
```
**other option** install with conda: https://anaconda.org/compbiocore/gsnap

2. Follow directions from http://research-pub.gene.com/gmap/ and http://research-pub.gene.com/gmap/src/README to install gsnap on Atlas with C. Moved gmap tool file to `/project/dot_files/software/` and also made link to home directory.
```
./configure --prefix=/location/to/gmap-2020-12-17 --with-gmapdb=/location/you/want/gmapsoftware
make
make check # run checks
make install
```
* Any C program with a Makefile can usually be installed with `./configure -> make -> make install`.

* If we look at the README instructions, it shows the configure -> make -> make install method of building. Here, we modified the `./configure` slightly to place the gmap executable in the software folder.

* executables including gmap are in `bin/`: `/project/dot_files/software/gmap-2020-12-17/bin`
</details>

3. See Jennifer's notes on local install of **featureCounts** and **samtools**. They were too difficult to install on Atlas, so will go with Miniconda instead.


## Setting up additional items on Atlas
1. Created several shortcuts to set up and run `debug` and `myjobs` via editing `.bashrc` file (rather than typing the long commands for each)
```
#debug to test code
salloc -N 1 -p atlas -t 01:00:00 --account=PROJnameHERE
#myjobs - see all jobs under UserName
squeue | grep UserNameHere
```

2. Run `source ~/.bashrc` to "restart" shell and apply latest changes to `.bashrc`. Atlas doesn't do this when you login, but Ceres does.

3. Created `inbox` and `outbox` so that we can do ssh mylocalfile.tar.gz UserName@atlas:inbox/. and don't have to think about 5gb home directory limit

## Restart: install everything in miniconda
Decided miniconda was the best approach to run gmap, featureCounts, samtools. Local install of these tools were too difficult (see Jennifer's notes).

1. Install miniconda to atlas by running `bash Miniconda3-latest-Linux-x86_64.sh`. Continuously press enter, even on prompt asking where to install miniconda3. You can move the source folder afterwards, in which case, it was moved to `project/`  and linked to `home/software/`.

2. Install samtools, gmap, subread (featureCounts) by first creating an environment file `gsnap_env.yml` in `miniconda3/envs/`
```
name: gsnap_env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.8
  - gmap
  - samtools
  - subread
```

3. From the command line (in `salloc` mode), create miniconda environment in same directory.
```
conda env create -f gsnap_env.yml
```

4. Check that environment was created:
```
conda env list          #<= list all environments
#> conda environments:
#> base                  *  /home/miniconda3
#> gsnap_env                /home/miniconda3/envs/gsnap_env
```

5. Activate conda environment and do a version check as a test that everything is working
```
conda activate gsnap_env
samtools --version        # check samtools version: 1.11
gmap --version        # gmap version: 2020-10-14
featureCounts -v    #featureCounts version: v2.0.1
```

6. To activate conda environment, activate local miniconda, and then `gsnap_env`.
```
#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=gsnap
#SBATCH --out=stdout.%j.%N.%x
#SBATCH --error=stderr.%j.%N.%x
#SBATCH --mail-user=myemail@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=ProjectName
set -e
set -u
set +eu
source /home/kathy.mou/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env
## gsnap/samtools/featureCount commands here
```

## Alignment with gsnap
### About gsnap
1. Workflow with gsnap: https://hbctraining.github.io/Intro-to-rnaseq-hpc-gt/lessons/08_rnaseq_workflow.html
* Set up directory tree like one below
```
rnaseq/
	├── raw_data/
	├── meta/
	├── results/
  |   |── gsnap/
	├── scripts/
	└── logs/
```

2. Literature about gsnap: https://link.springer.com/protocol/10.1007%2F978-1-4939-3578-9_15

3. The general format of the pipeline is:
```
(1) index genome -> (2) map reads to genome -> (3) get counts
```
* The counts will be sent to Diffential Expression analysis programs.

4. First step: Pre-process reference genome to create a genome index.
```
gmap_build -d <genome name> <path to genome fasta file>
```
* default value for -k is 15 (from https://github.com/juliangehring/GMAP-GSNAP/blob/master/README)

5. Second step: Map reads to genome
```
gsnap -d <genome> <read1_file> <read2_file>
```
* <genome> is the name of the genome database created by gmap_build

### Maize -- local gsnap
<details><summary>Notes of running **GSNAP** locally - made genomic index using local gsnap, continue the rest with miniconda3</summary>

1. Pre-process Zea mays B73 reference genome to create a genome index.
```
/project/projectdirectory/mydirectory/dot_files/software/gmap-2020-12-17/bin/gmap_build -d b73 /project/projectdirectory/mydirectory/rnaseq/maize/reference_genome/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
```
```
Submitted batch job 130262
```

2. Seems the command worked, but `b73/` ended up in `/project/projectdirectory/mydirectory/dot_files/software/gmapdb`. I looked at stderr file and noticed these two error messages found at the beginning and end of file
```
#lines 3-4
Destination directory not defined with -D flag, so writing files under /project/projectdirectory/mydirectory/dot_files/software/gmapdb
#line 338
/var/spool/slurmd/job130262/slurm_script: line 33: Wed: unbound variable
```
* Jennifer pointed out that when we did `./configure --something_gmapdb=/project/.../software` it set the `-D` to be in `software/gmapdb`
* `-d` is the genome name and `-D`  is where to store it
* in the future, also make sure to add a `-D desired/directory/` to command
* ignore error message on line 33 of stderr file.

### Output files
```
b73.chromosome	    b73.contig.iit     b73.ref061regiondb	b73.version
b73.chromosome.iit  b73.genomebits128  b73.ref153offsets64meta
b73.chrsubset	    b73.genomecomp     b73.ref153offsets64strm
b73.contig	    b73.maps	       b73.ref153positions
```

3. Map reads to genome
```
gsnap -d <genome> <read1_file> <read2_file>
```
Parameters to add
* -D genome directory
* -d genome database (is it `b73/`?)
* what else?

```
gsnap -d /project/projectdirectory/mydirectory/dot_files/software/gmapdb/b73/ -D /project/projectdirectory/mydirectory/rnaseq/maize/gsnap/ /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/*_1.fastq.gz /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/*_2.fastq.gz
```

Example command for gsnap:
```
gsnap -d grch38_chr1 -D /gstore/scratch/hpctrain/chr1_reference_gsnap \
-t 6 -M 2 -n 10 -N 1 \
--quality-protocol=sanger -w 200000 --pairmax-rna=200000 \
-E 1 -B 2 \
-A sam raw_data/Mov10_oe_1.subset.fq | \
samtools view -bS - | \
samtools sort - \
 > results/gsnap/Mov10_oe_1.Aligned.sortedByCoord.out.bam
```
</details>

### Maize -- miniconda3
1. Mapping RNA-seq reads to B73
  ```
  set -e
  set -u
  set +eu
  source /home/kathy.mou/miniconda3/etc/profile.d/conda.sh
  conda activate gsnap_env
  # ==== Mapping RNA-seq reads. Use miniconda3
  gsnap -d b73 -D /project/projectdirectory/mydirectory/dot_files/software/gmapdb/ \
  -t 6 -M 2 -n 10 -N 1 \
  --quality-protocol=sanger -w 200000 --pairmax-rna=200000 -E 1 -B 2 \
  -A sam /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/SRR1573504_1.fastq /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/SRR1573504_2.fastq| \
  samtools view -bS - | \
  samtools sort - \
  > /project/projectdirectory/mydirectory/rnaseq/maize/results/gsnap/SRR1573504_1_2.Aligned.sortedByCoord.out.bam
  ```
```
Submitted job 136482
```

2. `SRR1573504_1_2.Aligned.sortedByCoord.out.bam` file is empty so far...


### Bee -- miniconda3
1. Create genome index and map RNA-seq reads to B. impatiens
  ```
  set +eu
  source /home/mydirectory/miniconda3/etc/profile.d/conda.sh
  conda activate gsnap_env
  # ==== Connect the executable (either local or miniconda)
  # GMAP_BUILD=/project/projectdirectory/mydirectory/dot_files/software/gmap-2020-12-17/bin/gmap_build
  GMAP_BUILD=gmap_build
  # ==== Define input/output variables
  GENOME_NAME=B_impatiens
  GENOME_FASTA=/project/projectdirectory/mydirectory/rnaseq/bee/reference_genome_bee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa
  GMAPDB=/project/projectdirectory/mydirectory/dot_files/software/gmapdb
  # ==== Main Run
  ${GMAP_BUILD} -d ${GENOME_NAME} -D ${GMAPDB} ${GENOME_FASTA}
  # ==== Mapping RNA-seq reads. Use miniconda3
  gsnap -d ${GENOME_NAME} -D ${GMAPDB} \
  -t 6 -M 2 -n 10 -N 1 \
  --quality-protocol=sanger -w 200000 --pairmax-rna=200000 -E 1 -B 2 \
  -A sam /project/projectdirectory/mydirectory/rnaseq/bee/raw_data/1-A01-A1_S7_L002_R1_001.fastq | \
  samtools view -bS - | \
  samtools sort - \
  > /project/projectdirectory/mydirectory/rnaseq/bee/results/gsnap/1-A01-A1_S7_L002_R1.Aligned.sortedByCoord.out.bam
  ```
```
Submitted job 136483
```
Job seemed to complete. Resulted in `1-A01-A1_S7_L002_R1.Aligned.sortedByCoord.out.bam` in `/projectdirectory/mydirectory/rnaseq/bee/results/gsnap`

#### Output files
* `gmap_build`
```
b73.chromosome	    b73.contig.iit     b73.ref061regiondb	b73.version
b73.chromosome.iit  b73.genomebits128  b73.ref153offsets64meta
b73.chrsubset	    b73.genomecomp     b73.ref153offsets64strm
b73.contig	    b73.maps	       b73.ref153positions
```

* `gsnap`



## Differential expression with DESeq2
### Maize
### Bee
#### Output files
