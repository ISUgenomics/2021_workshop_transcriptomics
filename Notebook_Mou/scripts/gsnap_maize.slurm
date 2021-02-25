#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=b73
#SBATCH --out=stdout.%j.%N.%x
#SBATCH --error=stderr.%j.%N.%x
#SBATCH --mail-user=myem@il.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=ProjectName
set -e
set -u

# ==== Activate miniconda
set +eu
source /home/k/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env

# ==== Define input/output variables
REF_FASTA=/home/k/rnaseq/maize/reference_genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
REF_GFF=/home/k/rnaseq/maize/reference_genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff
GMAPDB=/project/f/k/dot_files/software/gmapdb
DB_NAME="b73" #<= was NAMV5, had forgotten to change this to b73

# === Set working directory and in/out variables
cd /project/f/k/rnaseq/maize/results/

# ==== Index Genome
#gmap_build -d ${DB_NAME} -D ${GMAPDB} ${REF_FASTA}

# ==== Mapping RNA-seq reads
for file in /home/k/rnaseq/maize/raw_data/*_1.fastq
  do
    file2=$(echo $file|sed -r 's/\_1/\_2/g')
    OUTFILE=$(basename ${file} | sed 's/_1\.fastq$//g')
    #echo ${OUTFILE}
    gsnap -d ${DB_NAME} -D ${GMAPDB} \
    -N 1 -t 16 -B 4 -m 5 \
    --input-buffer-size=1000000 \
    --output-buffer-size=1000000 \
    -A sam \
    --split-output=${DB_NAME}_${OUTFILE} ${file} ${file2} | \
    #echo $file
    #echo $file2
    samtools view --threads 16 -bS - > ${OUT_BAM}

# ==== Gene Counts
    OUT_COUNTS=${OUTFILE}_genecounts.txt
    featureCounts -T 16 -t gene -g ID \
      -a ${REF_GFF} \
      -o ${OUT_COUNTS} \
      ${OUT_BAM}
done
