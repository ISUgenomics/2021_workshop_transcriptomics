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
source /home/kathy.mou/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env


start=`date +%c`

# === Load Modules here
# module load <name of module>

# === Set working directory and in/out variables
cd /project/projectdirectory/mydirectory/rnaseq/maize/gsnap/

# === Main Program
# ==== Define Variables (this way you can swap out directories/input file)
#GSNP_DIR=/project/projectdirectory/mydirectory/dot_files/software/gmap-2020-12-17/bin
#REF_FILE=/project/projectdirectory/mydirectory/rnaseq/maize/reference_genome/Zea_mays.B73_RefGen_v4.dna.toplevel.fa

# ==== Run gmap_build with Maize Reference Genome B73 to make genomic index. Used natively installed gsnap
#${GSNP_DIR}/gmap_build -d b73 ${REF_FILE}

# ==== Mapping RNA-seq reads. Use miniconda3
gsnap -d b73 -D /project/projectdirectory/mydirectory/dot_files/software/gmapdb/ \
-t 6 -M 2 -n 10 -N 1 \
--quality-protocol=sanger -w 200000 --pairmax-rna=200000 -E 1 -B 2 \
-A sam /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/SRR1573504_1.fastq /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/SRR1573504_2.fastq| \
samtools view -bS - | \
samtools sort - \
> /project/projectdirectory/mydirectory/rnaseq/maize/results/gsnap/SRR1573504_1_2.Aligned.sortedByCoord.out.bam



end=`date +%c`
