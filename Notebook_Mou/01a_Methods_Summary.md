# Methods Summary

Written summary of methods performed in this repo.

## 1) Bee and Maize Data

### 1a) Table summary of bee and maize sequence data
| Sample | Reference Genome source | Notes |
| -- | -- | -- |
| Bumble bee (*Bombia impatiens*) | [Bombia impatiens](https://hymenoptera.elsiklab.missouri.edu/genome_fasta): Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa | Goal: Compare RNA transcriptome of bee brains living in normal vs heavy metal environments  <br/> Study details: Already have FastQC and multiQC reports of bee fastq files, QuantSeq 3'-generated data <br/> [Bumble bee metadata](/2021_workshop_transcriptomics/00a_Metadata.md) |
| [Maize Leaf](https://www.ebi.ac.uk/ena/browser/view/PRJNA260793) | [Maize reference (Zea mays B73)](ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.*gz) | Goal: "We utilized laser-microdissection RNAseq to identify genes that are differentially expressed along discrete cell/tissue-specific domains along the proximal-distal axis of wild-type leaf primordia undergoing ligule initiation, and compared transcript accumulation in wild type and liguleless1 mutant leaf primordia."" <br/> Study design: "In this study, we analysed the transcriptome associated with ligule formation using laser microdissection RNA-sequencing (LM-RNAseq). We quantified transcript accumulation in the PLB and adjacent pre-blade and pre-sheath regions of wild-type leaf primordia in order to identify candidate genes involved in proximal-distal patterning at the blade-sheath boundary. We also compared transcript accumulation in lg1-R mutants and wild-type siblings to identity genes acting downstream of LG1." <br/> Liguleless1 mutants lack ligules and auricles <br/> Reference: https://www.ebi.ac.uk/ena/browser/api/xml/SRP047035 <br/> [Maize metadata](/2021_workshop_transcriptomics/00a_Metadata.md) |

### 1b) Table summary of bee metadata
| group \ nest| A | B | C | D | E | F | total |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|Control| 18 | - | 8 | - | - |- |26|
|Exposed| - | 8 | - | 9 | 7 | 10 |34|

### 1c) Table summary of maize metadata
| group \ tissue | leaf blade (B) | leaf ligule (L) | leaf sheath (S) | all cell layers |
|:--|:-:|:-:|:-:|:-:|
|Control - all cell layers, B/L/S, 3 replicates (designated 3, 4, or 5) | 3 | 3 | 3 | - |
| Wild-type tissue from P6 leaf primordia, all cell layers (wtL) | - | - | - | 3 |
| Liguleless1 mutant tissue from P6 leaf primordia, all cell layers (lg1) | - | - | - | 3|
|adaxial L1 epidermal layer, B/L/S, 3 replicates (designated 1, 2, or 3) (*_L1, where * = B, L, or S)| 3 | 3 | 3 | - |

### 1d) Additional Notes
<details><summary>Bee Study Design Notes</summary>

#### Bee Study Design
* All females, foraging range 2km, exposed group of bees treated with heavy metals, soil high in lead and 3 other heavy metals associated with soils in urban areas. How does urbanization affect environment and these bees
* From group discussions asking which bee annotation to use (NCBI or Hymenoptera): Amy Toth recommends using Hymenoptera Base
```
  For social insect genomes, usually the "Official Gene Set" is the one that is on Hymenoptera Base (Elsik Lab).  NCBI has their own annotation for each genome as well.  The NCBI annotations tend to have fewer genes but are usually very high quality annotations.  I think most people in the field go with the OGS for a given species unless they are comparing across species, in which case they might use NCBI for consistency.
```
</details>


<details><summary>Other Notes</summary>

* HiSat2 (gene-level alignment), Kallisto (transcript-level alignment)
* Nests are confounding factors because there is no nest that spans both treatments. Can't tease out nest from treatments.
* Is there clustering by nests? See if there's nesting effect
* Ctrl_NestA, Ctrl_NestC, Exposed_NestB, etc.
* See if there's difference between Control nest A, nest C; Exposed Nest B, and others.
</details>

## 2) Software Setup

### 2a) Data transfer to HPC (Atlas dtn node)
* `ssh username@Atlas-dtn.hpc.msstate.edu`
https://www.hpc.msstate.edu/computing/atlas/

<details><summary>Bee</summary>

1. Working directory on Atlas: `/projectdirectory/mydirectory/rnaseq/bee`

2. Download all data from google drive to local computer in new folder called `bee/`

3. Copy fastq files from raw data directories to `bee/`:
```
find . -name *.fastq -exec cp '{}' "./raw/" ";"`
```

4. Copy `.sam` and `Bombus_impatiens_*` files from mapping directories to `mapping/`:
```
find . -name *.sam -exec cp '{}' "./mapping/" ";"`
find . -name Bombus_impatiens_* -exec cp '{}' "./mapping/" ";"
```
  * Bee reference genome in here: `Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa`

5. Tar raw and mapping directories together
```
tar -zcvf bee.tar.gz bee/
```

6. Uploaded `bee.tar.gz` to Atlas. Takes about 1.5h for 22GB.
</details>

<details><summary>Maize</summary>

1. Working directory on Atlas: `/projectdirectory/mydirectory/rnaseq/maize`

2. Fetch maize reference genome:
```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.*gz
```

3. Slurm script for fetching maize sequences to Atlas. See Jennifer's `atlas_maizedata.slurm` for template
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
</details>

### 2b) Set up project directory

<details><summary>Atlas Project Directory</summary>

```
ProjectDirectory/
  |_MyDirectory/
      |_bee/
          |_bee.tar.gz
          |_logs/
          |_mapping/
          |_meta/
          |_outbox/
          |_raw_data/
              |_test/
          |_reference_genome_bee/
          |_results/
              |_gsnap/
              |_multiqc/
          |_scripts/
      |_maize/
          |_gsnap/
              |_b73_reference_gsnap/
          |_logs/
          |_metadata/
          |_outbox/
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
</details>

### 2c) Install miniconda and RNAseq tools
<details><summary>Procedures for setting up miniconda on Atlas</summary>
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
source /h/k/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env
## gsnap/samtools/featureCount commands here
```
</details>

## 3) QC of fastq sequences with fastQC
<details><summary>Bee</summary>
Already done by Toth group, sequences look good.
</details>

<details><summary>Maize</summary>
1. Ran the following slurm script on Atlas
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Maizefastqc
#SBATCH --out=stdout.%j.%N.%x
#SBATCH --error=stderr.%j.%N.%x
#SBATCH --mail-user=em@il.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=f
# Set working directory
cd /h/k/rnaseq/maize/raw_data/
module load fastqc
fastqc -t 16 *.fastq
```

Moved fastqc files to subdirectory `fastqc/`

### Output files:
* `*.fastqc.html`
* `*.fastqc.zip`
</details>

## 4) Alignment with gsnap
<details><summary>About gsnap</summary>
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
</details>

<details><summary>Bee</summary>

### Bee -- miniconda3
1. Create genome index and map RNA-seq reads to *B. impatiens*
  ```
  set +eu
  source /h/mydirectory/miniconda3/etc/profile.d/conda.sh
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

### Output file:
* `1-A01-A1_S7_L002_R1.Aligned.sortedByCoord.out.bam` in `/projectdirectory/mydirectory/rnaseq/bee/results/gsnap`

2. Run featureCounts
```
#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=bee
#SBATCH --out=stdout.%j.%N.%x
#SBATCH --error=stderr.%j.%N.%x
#SBATCH --mail-user=myem@il.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=fsepru
set -e
set -u
# ==== Activate miniconda
set +eu
source /h/k/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env
# ==== Define input/output variables
REF_NAME=B_impatiens
#REF_FASTA=/project/f/k/rnaseq/bee/reference_genome_bee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa
GMAPDB=/project/f/k/dot_files/software/gmapdb
REF_GFF=/project/f/k/rnaseq/bee/reference_genome_bee/GCF_000188095.3_BIMP_2.2_genomic.gff.gz
# === Set working directory and in/out variables
cd /project/f/k/rnaseq/bee/results/
# ==== Gene Counts
for FILE in /h/k/rnaseq/bee/raw_data/test/*.bam
  do
    OUT_COUNTS=${FILE}_genecounts.txt
    OUTBAM=${FILE}

    featureCounts -T 16 -t gene -g ID \
       -a ${REF_GFF} \
       -o ${OUT_COUNTS} \
       ${OUTBAM}
done
```

### Output files
* `gmap_build`
```
B_impatiens.chromosome	    B_impatiens.contig.iit     B_impatiens.ref061regiondb	B_impatiens.version
B_impatiens.chromosome.iit  B_impatiens.genomebits128  B_impatiens.ref153offsets64meta
B_impatiens.chrsubset	    B_impatiens.genomecomp     B_impatiens.ref153offsets64strm
B_impatiens.contig	    B_impatiens.maps	       B_impatiens.ref153positions
```

* `gsnap`

```
*.fastq.Aligned.sortedByCoord.out.bam
*.fastq
```

* `featureCounts`

```
*genecounts.txt
*genecounts.txt.summary
```
</details>

<details><summary>Maize</summary>

### Maize -- miniconda3
1. Mapping RNA-seq reads to B73
  ```
  set -e
  set -u
  set +eu
  source /h/k/miniconda3/etc/profile.d/conda.sh
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

2. Run gsnap (adapted from `2021_workshop_transcriptomics/Notebook_Severin/Maize/02_gsnap.md` `gsnapScript.sh` and Jennifer's `Maize_Runner.slurm`):
```
#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Maize
#SBATCH --out=stdout.%j.%N.%x
#SBATCH --error=stderr.%j.%N.%x
#SBATCH --mail-user=myem@il.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=f
set -e
set -u
start=`date +%s`
# === Load Modules here and link executables
# = Atlas HPC
set +eu
source /h/k/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env
GMAP_BUILD=gmap_build
GSNAP=gsnap
SAMTOOLS=samtools
FEATURECOUNTS=featureCounts
# === Set working directory and in/out variables
cd /project/f/k/rnaseq/maize/results/
# === Input / Output Variables
REF_NAME=b73
REF_FILE=/h/k/rnaseq/maize/reference_genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
REF_GFF=/h/k/rnaseq/maize/reference_genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff
GMAPDB=/project/f/k/dot_files/software/gmapdb
# # === Main Program
# (1) Index Genome
#${GMAP_BUILD} \
#  --gunzip \
#  -d ${REF_NAME} \
#  -D ${GMAPDB} \
#  ${REF_FILE}
for FILE in /h/k/rnaseq/maize/raw_data/*_1.fastq
do
  READ_NAME=$(basename ${FILE} | sed 's:_1.fastq::g')
  DIR_NAME=$(dirname ${FILE})
  READ_R1=${DIR_NAME}/${READ_NAME}_1.fastq
  READ_R2=${DIR_NAME}/${READ_NAME}_2.fastq
  OUT_BAM=${READ_NAME}.aligned.out.bam
  OUT_COUNTS=${READ_NAME}_genecounts.txt
  echo "Processing ... ${READ_NAME}"
# (2) Map Reads:
  ${GSNAP} \
    --gunzip \
    -d ${REF_NAME} \
    -D ${GMAPDB} \
    -N 1 -t 16 -B 4 -m 5 \
    --input-buffer-size=1000000 \
    --output-buffer-size=1000000 \
    -A sam \
    ${READ_R1} ${READ_R2} | \
    ${SAMTOOLS} view --threads 16 -bS - > ${OUT_BAM}
# (3) Get feature counts
  ${FEATURECOUNTS} -T 16 -t gene -g ID \
    -a ${REF_GFF} \
    -o ${OUT_COUNTS} \
    ${OUT_BAM}
done
end=`date +%s`
# === Log msgs and resource use
scontrol show job ${SLURM_JOB_ID}
echo "ran Bee_Runner.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
```

### Output files:
* `*genecounts.txt`
* `*genecounts.txt.summary`
</details>

## 5) Counts with featureCounts

<details><summary>Useful references for understanding how FeatureCounts works</summary>

* FeatureCounts User Guide: http://www.bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
* Nice diagram showing cases of the effect of `countmultioverlap` on `overlapmethod`: https://www.mathworks.com/help/bioinfo/ref/featurecount.html
* Reads that map to multiple transcripts or align ambigiuously and not sure what strand they are (+ / -) aka ambiguous - problem with short reads
* can ignore multimapped duplicates, ignore duplicates - unique reads will be based on read count.
* 7th column is the read count
* Apparently you can run featureCounts on sam files (tested by Andrew Severin)
* `combine.R` notes:
  * `(featureCount_files <- list.files(path = dir_org, pattern = "*genecounts.txt$", full.names = TRUE))`
    * This includes "/" and directory of *genecounts.txt file. $ = end of string (regular expression). The outside () prints out object (`featureCount_files`)
  * `writexl::write_xlsx` = package::function, writexl is the package, write_xlsx is the function within writexl package. Simple way to call a function without having to hassle with installing and loading new packages.
</details>

<details><summary>Bee</summary>

1. Downloaded featureCount output files from Atlas.

2. Run featureCounts output through `combine.R`.

```
#! /usr/bin/env Rscript
# Auth: Jennifer Chang & Mou
# Date: 2021/03/04
# Desc: Combine featureCounts output (1st and last column) files for Bee and Maize. The text files (*.genecounts.txt) generated from this script will be used for DESeq2.

# === Load Libraries
library(tidyverse)
library(magrittr)
library(readxl)

###### BEE #######

# === Get list of featureCount output files
dir_org="~/Desktop/bee/"        # counts are in a "bee" or "maize" subdirectory
(featureCount_files <- list.files(path = dir_org, pattern = "*genecounts.txt$", full.names = TRUE)) #includes "/" and directory of *genecounts.txt file, $ = end of string
#(featureCount_files2 <- list.files(path = dir_org, pattern = "*genecounts.txt$")) #only file name

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" )  %>%
  select(Geneid, ends_with(".bam"))  %>%              # Get 1st and last column (column was named after bam file)
  pivot_longer(cols=ends_with(".bam")) %>%           # Melt data (tidy data)
  mutate(
    name = gsub(".Aligned.sortedByCoord.out.bam", "", name)        # No longer need the bam extension, easier to read
  ) %>%
  mutate(name = gsub("/home/kathy.mou/rnaseq/bee/raw_data/test/", "", name)        # No longer need the bam extension, easier to read
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(Geneid, ends_with(".bam")) %>%
    pivot_longer(cols=ends_with(".bam")) %>%
    mutate(
      name = gsub(".Aligned.sortedByCoord.out.bam", "", name)
    ) %>%
    mutate(name = gsub("/home/kathy.mou/rnaseq/bee/raw_data/test/", "", name))       
  data = rbind(data, temp)
}

# === Convert to excel like data (wider)
wide_data <- data %>%
  pivot_wider(id_cols=Geneid)

# === Save tab delimited file (smaller file size)
write_delim(wide_data,
            paste(dir_org, "/bee.genecounts.out.txt", sep = ""),
            delim="\t")

# === Save Excel file (can be easier to work with)
writexl::write_xlsx(wide_data,
                    path=paste(dir_org, "/bee.genecounts.xlsx", sep = ""))
#package :: function
```

### Output files:
* bee.genecounts.out.txt
* bee.genecounts.xlsx
</details>

<details><summary>Maize</summary>

1. Downloaded featureCount output files from Atlas.

2. Run featureCounts output through `combine.R`.

```
#! /usr/bin/env Rscript
# Auth: Jennifer Chang & Mou
# Date: 2021/03/04
# Desc: Combine featureCounts output (1st and last column) files for Bee and Maize. The text files (*.genecounts.txt) generated from this script will be used for DESeq2.

# === Load Libraries
library(tidyverse)
library(magrittr)
library(readxl)

###### MAIZE #######

# === Get list of featureCount output files
dir_org="~/Desktop/maize/"        # counts are in a "bee" or "maize" subdirectory
(featureCount_files <- list.files(path = dir_org, pattern = "*genecounts.txt$", full.names = TRUE)) #includes "/" and directory of *genecounts.txt file, $ = end of string
#(featureCount_files2 <- list.files(path = dir_org, pattern = "*genecounts.txt$")) #only file name

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" )  %>%
  select(Geneid, ends_with(".bam"))  %>%              # Get 1st and last column (column was named after bam file)
  pivot_longer(cols=ends_with(".bam")) %>%           # Melt data (tidy data)
  mutate(
    name = gsub(".aligned.out.bam", "", name)        # No longer need the bam extension, easier to read
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(Geneid, ends_with(".bam")) %>%
    pivot_longer(cols=ends_with(".bam")) %>%
    mutate(
      name = gsub(".aligned.out.bam", "", name))       
  data = rbind(data, temp)
}

# === Convert to excel like data (wider)
wide_data <- data %>%
  pivot_wider(id_cols=Geneid)

# === Save tab delimited file (smaller file size)
write_delim(wide_data,
            paste(dir_org, "/maize.genecounts.out.txt", sep = ""),
            delim="\t")

# === Save Excel file (can be easier to work with)
writexl::write_xlsx(wide_data,
                    path=paste(dir_org, "/maize.genecounts.xlsx", sep = ""))
```

### Output files:
* maize.genecounts.out.txt
* maize.genecounts.xlsx
</details>

## 6) QC of featureCounts summary files with MultiQC
<details><summary>About featureCounts and MultiQC</summary>

* Check out this [link](https://multiqc.info/docs/) on how to run multiQC with featureCounts `*.summary` files. This is another way to assess read alignment quality.
* Read MultiQC to assess read alignment (can try this on gsnap output): http://www.bea.ki.se/documents/Intro2RNAseq.pdf
* What parameter to use in gsnap to remove bad reads? (unassigned: unmapped, ambiguous, multimapping, no features) - don't need to, if run default, the uniquely read counts will be the mapped reads.
* Make sure to set package cache in desired project directory. Go to `.condarc` file in home directory and  modify to something like this:
  ```
  pkgs_dirs:
  - /p/f/k/my_pkg_cache
  ```
</details>

<details><summary>Bee</summary>

1. Make `/p/f/k/rnaseq/bee/results/multiqc/` directory and copy `/p/f/k/rnaseq/bee/raw_data/testing/*.summary` files to `multiqc/`. Run `multiqc .`

2. Check out `bee.multiqc_report.html`.
Notice that `1-E07-F5_S61_L002_R1_001` had the lowest number of assigned reads (223,397). All others had at least 1M reads.

3. Use MultiQC Toolbox on html page to export `featureCounts_assignment_plot` image and save as `Bee_featureCounts_multiqc_plot.png`
![](results/Bee_featureCounts_multiqc_plot.png)<!-- -->
### Output files:
* `multiqc_data/` <= for both bee and maize
* `bee.multiqc_report.html`
</details>

<details><summary>Maize</summary>

1. Make `/p/f/k/rnaseq/maize/results/multiqc/` and copy `/p/f/k/rnaseq/maize/results/*.summary` files to `multiqc/`. Run `multiqc .`

2. Check out  `maize.multiqc_report.html`. Note that  SRR1573520 has a lot of unassigned_multimapping reads.

3. Use MultiQC Toolbox on html page to export `featureCounts_assignment_plot` image and save as `maize_featureCounts_multiqc_plot.png`
![](results/maize_featureCounts_multiqc_plot.png)<!-- -->
### Output files:
* `multiqc_data/`
* `maize.multiqc_report.html`
</details>

## 7) Differential expression with DESeq2
<details><summary>Notes about DESeq2</summary>

* Sathesh says basemeans correlate with read counts: larger basemean values = more read counts
* DESeq2 reference guide: http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink
</details>

<details><summary>Bee</summary>
1. Run `bee_maize_deseq2.Rmd`.

2. Generate `Bee_AllExposedvsAllControlGene.csv` file, saved in `results/`. Made notes in Rmd file.

3. A few things to look for:
  * `eval=FALSE` was added within `{r}` in each code chunk so that when I ran Knit, I don't have to have maize data to generate markdown file and output.
  * line 71: try boxplot with group (GSNAP RNA Gene Profiles) as "x" within ggplot function
  * take out 1_E07 sample since it had very low reads??
  * line 141: test out the mutate function when creating `meta_df` object
  * line 167: test out `res <- res[order(res$padj), ]`
  * line 183:
  ```
  PCA by Nest + treatment
ANOVA stats to look at treatment, nest effects on variation (only genes with large fold-change, p<0.05)
  ```
  * line 207: expression = normalized read count? Why does 1-B11 have super high #? Same with 1-E07
3. Ran a couple diagnostic plots: linechart and violin plots of count data
![](results/Bee_GSNAP_RNAseq_CountProfiles_Linechart.png)<!-- -->

![](results/Bee_GSNAP_RNAseq_LogCountProfiles_ViolinPlot.png)<!-- -->

4. Made PCA of all genes , subset by Treatment group (Exposed vs Control). Looks exactly like Rick's PCA except flipped upside down.
![](results/Bee_AllGenes_PCA_ExposedvsControl.png)<!-- -->

5. Subset by nest `Bee_Nest_DESeq2.csv`, but list of genes and p-values do not appear to be different from `Bee_AllExposedvsAllControlGene.csv`.

6. Made PCA of all genes, subset by nest. Didn't see any particular clustering by nest.
![](results/Bee_AllGenes_PCA_Nest.png)<!-- -->

7. Made volcano plot of all genes to see a big picture of proportion of not differentially regulated genes vs up- or down-regulated genes
![](results/Bee_AllGenes_VolcanoPlot.png)<!-- -->

8. Made heatmap of bee count matrix
![](results/Bee_HeatmapOfCountMatrix.png)<!-- -->

9. Made heatmap of bee sample-to-sample distances
![](results/Bee_HeatmapOfSampleToSampleDistances.png)<!-- -->

10. Things to try
* Apply count outlier detection with Cook's distance `res$stat`

### Output files
* `Bee_AllExposedvsAllControlGene.csv`
* `Bee_Nest_DESeq2.csv`
</details>

<details><summary>Maize</summary>

1. Ran `bee_maize_deseq2.Rmd` up to adding metadata csv file. Found maize metadata here: https://www.ebi.ac.uk/ena/browser/view/PRJNA260793. Downloaded report (tsv file). Most important columns are run_accession (sample IDs) and sample_title (groups). Saved as `maize_metadata_All_Info.csv`.
2. Made a new metadata file `maize_metadata.csv` so that it only includes the columns `run_accession` and `Tissue`. `Tissue` derived from `sample_title`. I converted like the following:
| sample_title | Tissue |
| -- | -- |
| B-3 | B |
| L-3 | L |
| S-3 | S|
| wtL-1 | wtL |
| lg1-1 | lg |
| B_L1.1 | B_L |
| L_L1.1 | L_L |
| S_L1.1 | S_L |
3. Ran a couple diagnostic plots: linechart and violin plots of count data
![](results/Maize_GSNAP_RNAseq_CountProfiles_Linechart.png)<!-- -->

![](results/Maize_GSNAP_RNAseq_LogCountProfiles_ViolinPlot.png)<!-- -->

4. Subset by tissue `Maize_Tissue_DeSeq2.csv`.

5. Made PCA of all genes, subset by Tissue type.
![](results/Maize_AllGenes_PCA_Tissue.png)<!-- -->

6. Make volcano plot of all genes to see a big picture of proportion of not differentially regulated genes vs up- or down-regulated genes


### Output files
* `Maize_Tissue_DeSeq2.csv`
</details>

## 8) Network analysis with WGCNA
* which count table to use?
* use list of differentially expressed genes only
* Gene ontology enrichment analysis: Pretty easy as long as you have gene annotation (gff?); blastp to most related organism (find most recent bumblebee ontology papers for GO assignments, blast results) , blast2GO
* other options: goseq (R package), find GO file from other bee RNAseq papers (are they pulling the most recent GO - go to their source code and modify if necessary)
* create a correlation network (spearman)
* Cytoscape
