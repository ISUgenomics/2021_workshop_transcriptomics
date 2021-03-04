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
