#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Bees
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=jenchang@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=isu_gif_vrsc

set -e
set -u

start=`date +%s`

# === Load Modules here and link executables

# = Nova HPC
# module load gmap-gsnap
# module load samtools
# module load subread

# = Atlas HPC
set +eu
source /home/jennifer.chang/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env
GMAP_BUILD=gmap_build
GSNAP=gsnap
SAMTOOLS=samtools
FEATURECOUNTS=featureCounts

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Input / Output Variables
REF_NAME=Bombus
REF_FILE=data_bee/ref/GCA_000188095.3_BIMP_2.1_genomic.fna.gz
REF_GFF=data_bee/ref/GCA_000188095.3_BIMP_2.1_genomic.gff.gz
GMAPDB=gmapdb
# See forloop for the directory of reads

# === Main Program
# (1) Index Genome
${GMAP_BUILD} \
  --gunzip \
  -d ${REF_NAME} \
  -D ${GMAPDB} \
  ${REF_FILE}

for FILE in data_bee/reads/*.fastq
do
  READ_NAME=$(basename ${FILE} | sed 's:_L002_R1_001.fastq::g')
  DIR_NAME=$(dirname ${FILE})
  READ_R1=${DIR_NAME}/${READ_NAME}_L002_R1_001.fastq
  OUT_BAM=${READ_NAME}.aligned.out.bam
  OUT_COUNTS=${READ_NAME}_genecounts.txt

# (2) Map Reads
  ${GSNAP} \
    -d ${REF_NAME} \
    -D ${GMAPDB} \
    -N 1 -t 16 -B 4 -m 5 \
    --input-buffer-size=1000000 \
    --output-buffer-size=1000000 \
    -A sam \
    ${READ_R1} | \
    ${SAMTOOLS} view --threads 16 -bS - > ${OUT_BAM}

# (3) Gene Counts
  ${FEATURECOUNTS} -T 16 -t gene -g ID \
    -a ${REF_GFF} \
    -o ${OUT_COUNTS} \
    ${OUT_BAM} 

done

end=`date +%s`

# === Log msgs and resource use                          
scontrol show job ${SLURM_JOB_ID}
echo "ran gsnap.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
