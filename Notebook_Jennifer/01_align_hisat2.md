# Hisat2 - alignment

Authors: Siva & Jennifer

Log onto Nova HPC and check if the following programs are installed

* <b>fastqc</b> - to quality check the fastq files
* <b>multiqc</b> - to merge the fastqc reports into one html
* <b>HISAT2</b> - for alignment
* <b>parallel</b> - to run some of the above commands in parallel

## Check Nova HPC Modules

The usual way to check if a program is installed is using:

```
(module avail) &> module_list.txt        # <= List all modules
grep -i hisat module_list.txt            # <= pull out any modules with hisat in the name
```

The following will be in our slurm scripts:

```
module load hisat
```

<details><summary>Install HiSat2 - for Atlas later, not done</summary>

Visit the HISAT2 download page, fetch the linux version.

* [http://daehwankimlab.github.io/hisat2/download/](http://daehwankimlab.github.io/hisat2/download/)

Maybe secure copy (`scp`) it over to Atlas HPC...

```
scp ~/Downloads/hisat2-2.2.1-Linux_x86_64.zip atlas:inbox/.
```

Heh, maybe I should add a discussion of `~/.ssh/config` files. It shortens the HPC address, easier to navigate to.

</details>

## Datasets

Andrew and Siva loaded the datafiles onto Nova in the following location.

```
/work/gif/TranscriptomicsWorkshop/
   |_ bumbleBee/
      |_ 01_Genome/
      |  |_ GCF_000188095.3_BIMP_2.2_feature_count.txt.gz
      |  |_ GCF_000188095.3_BIMP_2.2_genomic.fna
      |  |_ GCF_000188095.3_BIMP_2.2_genomic.gff.gz
      |
      |_ 1-A01-A1_S7_L002_R1_001.fastq.gz
      |_ ... rest of reads
      
   |_ Maize/
      |_ 01_Genome/
      |   |_ GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
      |   |_ GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff
      |
      |_ SRR1573504_1.fastq
      |_ SRR1573504_2.fastq
      |_ ... rest of paired reads
```

# 01) Quality Check

Following the commandßs from Siva's Notebook, run `fastqc` and `multiqc` to quality check the input fastq files.

```
# === Modules
module load fastqc/0.11.7-3flwcvl
module load py-multiqc/1.5-py2-lqqx3ht
#module load singularity/3.5.2-py2-ewuuq5o

# === Input / Output
INDIR=Maize
OUTDIR=01_QC

# === Main Program
parallel -j32 "fastqc -o ${OUTDIR}" ::: ${INDIR}/*.fastq
cd 01_QC
multiqc

#singularity exec --bind $PWD /work/gif/Siva/multiqc_latest.sif multiqc.
```



# 1) Hisat2 Index the reference genomes

Following the commnds from Siva's Notebook. 

```
# === Modules
module load hisat2/2.2.0-5kvb7f2

# === Input / Output variables
REF_FILE=Zm_B73-REFERENCE_Nam.5.0.fa
REF_NAME=b73
LOGFILE=v5_build.log

# === Main Program

hisat2-build -p 18 \
  ${REF_FILE} \
  ${REF_NAME} \
  >& ${LOGFILE}
```

# 2) Hisat2 Map the reads

Following the commands from Siva's Notebook.

```
hisat2 -p 8 \
  -x ${REF_NAME} \
  -1 ${READ_R1} \
  -2 ${READ_R2} \
  -S ${READ_NAME}.sam \
  2> ${READ_NAME}.log
  
samtools view \
  --threads 8 \
  -bS \
  -o ${READ_NAME}.bam

```

Following commands in Siva's Notebook.

HiSat2 alignment followed by Stringtie/Ballgown. This is definitely not done! Writing up notes before I start running commands.


## Maize

Hmm, does HISAT have a maize or bee genome? I should be able to build it manually right?

1. Index Genome

```
# === Input / Output Variables
REF_FILE=data_maize/ref/*.fna.gz
REF_NAME=b73

# === Main Program
hisat2-build ${REF_FILE} ${REF_NAME}
```

Will generate genome index files with a `*.ht21` file extension.

2. Run HiSAT to get counts?

Hmm... not a fan of this input format... maybe there's a param I'm missing.

Still reading through the documentation...

* [http://daehwankimlab.github.io/hisat2/manual/](http://daehwankimlab.github.io/hisat2/manual/)

```
# === Input / Output Variables
REF_FILE=data_maize/ref/*.fna.gz
REF_NAME=b73

# === Main Program
# create comma separated list of left and right reads...
ls data_maize/reads/*_1.fastq.gz |tr '\n' ',' > reads_1.txt
ls data_maize/reads/*_2.fastq.gz |tr '\n' ',' > reads_2.txt

hisat2 -x ${REF_NAME} \
  -1 $(cat reads_1.txt) \
  -2 $(cat reads_2.txt) \
  -p 16 \
  -S ${READ_NAME}.sam
```

Can I redirect this to `samtools view` to get a smaller bam file? Wait is this combining all together? Maybe I still need the file loop.
