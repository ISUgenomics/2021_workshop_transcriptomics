#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --mem 300000
#SBATCH --job-name=MAP
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=jenchang@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
# --account=project_name

set -e
set -u

start=`date +%s`

# === Load Modules here
# module load <name of module>
module load gmap-gsnap
module load samtools

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Input / Output Variables
REF_FILE=/work/gif/TranscriptomicsWorkshop/jennifer/temp_data/Maize/ref/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
REF_NAME=b73
GMAPDB=/work/gif/TranscriptomicsWorkshop/jennifer/temp_data/gmapdb
READ_R1=/work/gif/TranscriptomicsWorkshop/jennifer/temp_data/Maize/reads/SRR1573504_1.fastq.gz
READ_R2=/work/gif/TranscriptomicsWorkshop/jennifer/temp_data/Maize/reads/SRR1573504_2.fastq.gz
READ_NAME=SRR1573504
OUT_BAM=${READ_NAME}.aligned.sortedByCoord.out.bam
# === Main Program

# Index Reference
# gmap_build -g -d ${REF_NAME} -D ${GMAPDB} ${REF_FILE}

# Map Reads (test one read)
gsnap -d ${REF_NAME} -D ${GMAPDB} --gunzip \
  -t 16 -M 2 -n 10 -N 1 \
  --quality-protocol=sanger -w 200000 --pairmax-rna=200000 \
  -E 1 -B 2 \
  -A sam \
  ${READ_R1} ${READ_R2} | \
  samtools view --threads 16 -bS - | \
  samtools sort -m 30G --threads 16 - > ${OUT_BAM}

end=`date +%s`

# === Log msgs and resource use                          
scontrol show job ${SLURM_JOB_ID}
echo "ran gsnap.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
