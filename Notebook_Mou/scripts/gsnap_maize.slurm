#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=gsnap
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

# ==== Connect the executable (either local or miniconda)
# GMAP_BUILD=/project/f/k/dot_files/software/gmap-2020-12-17/bin/gmap_build
GMAP_BUILD=gmap_build

# ==== Define input/output variables
GENOME_NAME=B_impatiens
GENOME_FASTA=/project/f/k/rnaseq/bee/reference_genome_bee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa
GMAPDB=/project/f/k/dot_files/software/gmapdb

# ==== Mapping RNA-seq reads. Use miniconda3
export GMAPDB=/project/f/k/dot_files/software/gmapdb/b73/
DB_NAME="NAMV5" #<= forgot to change this to b73
# Note: "-N" option for detecting novel splice sites, remove if not needed (0=OFF; 1=ON)
for file in /home/k/rnaseq/maize/raw_data/*_1.fastq
  do
    file2=$(echo $file|sed -r 's/\_1/\_2/g')
    OUTFILE=$(basename ${file} | sed 's/_1\.fastq$//g')
#echo ${OUTFILE}
    gsnap -d ${DB_NAME} -N 1 -t 8 -B 4 -m 5 --input-buffer-size=1000000 --output-buffer-size=1000000 -A sam --split-output=${DB_NAME}_${OUTFILE} ${file} ${file2}
#echo $file
#echo $file2
done
