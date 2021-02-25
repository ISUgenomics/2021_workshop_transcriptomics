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
source /home/k/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env

# ==== Define input/output variables
REF_NAME=B_impatiens
#REF_FASTA=/project/f/k/rnaseq/bee/reference_genome_bee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa
GMAPDB=/project/f/k/dot_files/software/gmapdb
REF_GFF=/project/f/k/rnaseq/bee/reference_genome_bee/GCF_000188095.3_BIMP_2.2_genomic.gff.gz
#OUTBAM=${READNAME}.Aligned.sortedByCoord.out.bam

# === Set working directory and in/out variables
cd /project/f/k/rnaseq/bee/results/

# ==== Index Genome
#${GMAP_BUILD} -d ${REF_NAME} -D ${GMAPDB} ${REF_FASTA}

# ==== Mapping RNA-seq Reads
#for readname in /home/k/rnaseq/bee/raw_data/test/*.fastq
#        do
#                OUTBAM=${readname}.Aligned.sortedByCoord.out.bam
#                gsnap -d ${REF_NAME} -D ${GMAPDB} \
#                        -t 6 -M 2 -n 10 -N 1 \
#                        --quality-protocol=sanger -w 200000 --pairmax-rna=200000 \
#                        -E 1 -B 2 \
#                        -A sam \
#                        ${readname} | \
#                        samtools view -bS - > ${OUTBAM}
#               echo ${readname}
#               echo ${OUTBAM}
#        done

# ==== Gene Counts
for FILE in /home/k/rnaseq/bee/raw_data/test/*.bam
  do
    OUT_COUNTS=${FILE}_genecounts.txt
    OUTBAM=${FILE}

    featureCounts -T 16 -t gene -g ID \
       -a ${REF_GFF} \
       -o ${OUT_COUNTS} \
       ${OUTBAM}
done
