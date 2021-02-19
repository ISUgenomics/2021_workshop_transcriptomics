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
#SBATCH --account=fsepru
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
#OUTBAM=${READNAME}.Aligned.sortedByCoord.out.bam

# ==== Main Run
${GMAP_BUILD} -d ${GENOME_NAME} -D ${GMAPDB} ${GENOME_FASTA}

# ==== Mapping RNA-seq reads. Use miniconda3
for readname in /home/k/rnaseq/bee/raw_data/test/*.fastq
        do
                OUTBAM=${readname}.Aligned.sortedByCoord.out.bam
                gsnap -d ${GENOME_NAME} -D ${GMAPDB} \
                        -t 6 -M 2 -n 10 -N 1 \
                        --quality-protocol=sanger -w 200000 --pairmax-rna=200000 \
                        -E 1 -B 2 \
                        -A sam \
                        ${readname} | \
                        samtools view -bS - > ${OUTBAM}
#               echo ${readname}
#               echo ${OUTBAM}
        done
