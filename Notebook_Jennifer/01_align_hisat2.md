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
module load samtools/1.10-py3-xuj7ylj
module load subread/1.6.0-ak6vxhs

# === Input / Output variables
REF_NAME=b73
REF_FILE=GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
REF_GFF=GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff
LOGFILE=${REF_NAME}_build.log

# processes or threads
PROC=16

# === Main Program

hisat2-build -p ${PROC} \
  ${REF_FILE} \
  ${REF_NAME} \
  >& ${LOGFILE}
```

Will generate genome index files with a `*.ht21` file extension.

# 2) Hisat2 Map the reads

Following the commands from Siva's Notebook. Can either use `parallel` (Siva) or `for FILE in folder/*.fastq` (Mou). 

```
PROC=16

for FILE in Maize/*_1.fastq
do
  READ_NAME=$(basename ${FILE} | sed 's:_1.fastq::g')
  READ_R1=${FILE}
  READ_R2=${READ_NAME}_2.fastq
  OUT_BAM=${READ_NAME}.aligned.out.bam

  # (2) Map reads to indexed genome
  hisat2 -p ${PROC} \
    -x ${REF_NAME} \
    -1 ${READ_R1} \
    -2 ${READ_R2} |\
    samtools view \
    --threads ${PROC} \
    -bS \
    -o ${OUT_BAM}
done
```

# 3) featureCounts, get gene counts

Same instructions as GSNAP

```
for FILE in *.bam
do
  OUT_COUNTS=$(basename ${FILE} | sed 's:.bam::g')_genecounts.txt
  
  featureCounts \
    -T ${PROC} \
    -t gene \
    -g ID \
    -a ${REF_GFF} \
    -o ${OUT_COUNTS}
done

```

HiSat2 alignment followed by Stringtie/Ballgown. This is definitely not done! Writing up notes before I start running commands.


# Maize - Hisat2 Run on Nova

<b>Maize_Runner_hisat2.slurm</b>

```
#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Index
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=jenchang@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
# --account=isu_gif_vrsc

set -e
set -u

start=`date +%s`

# === Load Modules here and link executables

# = Nova HPC
module load hisat2/2.2.0-5kvb7f2
module load samtools/1.10-py3-xuj7ylj
module load subread/1.6.0-ak6vxhs

# = Atlas HPC
# 

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Input / Output Variables
REF_NAME=b73
REF_FILE=Maize/01_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
REF_GFF=Maize/01_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff
BAMDIR=02_BAM_Maize
COUNTDIR=03_Counts_Maize
LOGFILE=${REF_NAME}.log
[[ -d ${BAMDIR} ]]   || mkdir -p ${BAMDIR}
[[ -d ${COUNTDIR} ]] || mkdir -p ${COUNTDIR}

HISAT2=hisat2
HISAT2_BUILD=hisat2-build
SAMTOOLS=samtools
FEATURECOUNTS=featureCounts

# === Main Program

# (1) Index Genome
${HISAT2_BUILD} \
  -p 16 \
  ${REF_FILE} \
  ${REF_NAME} \
  >& ${LOGFILE}

# === Switch this to the parallel command from Siva's notes
# 
for FILE in Maize/*_1.fastq
do

  READ_NAME=$(basename ${FILE} | sed 's:_1.fastq::g')
  DIR_NAME=$(dirname ${FILE})
  READ_R1=${DIR_NAME}/${READ_NAME}_1.fastq
  READ_R2=${DIR_NAME}/${READ_NAME}_2.fastq
  OUT_SAM=${READ_NAME}.aligned.out.sam
  OUT_BAM=${BAMDIR}/${READ_NAME}.aligned.out.bam
  OUT_COUNTS=${COUNTDIR}/${READ_NAME}_genecounts.txt
  echo "Processing ... ${READ_NAME}"

# (2) Map Reads:
  ${HISAT2} \
    -p 16 \
    -x ${REF_NAME} \
    -1 ${READ_R1} \
    -2 ${READ_R2} |\
    ${SAMTOOLS} view --threads 16 -bS - > ${BAMDIR}/${OUT_BAM}

#    -S ${OUT_SAM}  | \
# (3) Get feature counts
  ${FEATURECOUNTS} -T 16 -t gene -g ID \
    -a ${REF_GFF} \
    -o ${OUT_COUNTS} \
    ${OUT_BAM} 

done

end=`date +%s`

# === Log msgs and resource use                          
scontrol show job ${SLURM_JOB_ID}
echo "ran Maize_Runner_hisat2.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
```

# Bee - Hisat2 Run on Nova

Either link or copy and paste here.

```
```
