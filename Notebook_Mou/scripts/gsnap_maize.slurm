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

start=`date +%c`

# === Load Modules here
# module load <name of module>

# === Set working directory and in/out variables
cd /project/projectdirectory/mydirectory/rnaseq/maize/gsnap/

# === Main Program

# Maize Reference Genome B73: make genomic index
/project/projectdirectory/mydirectory/dot_files/software/gmap-2020-12-17/bin/gmap_build -d b73 /project/projectdirectory/mydirectory/rnaseq/maize/reference_genome/Zea_mays.B73_RefGen_v4.dna.toplevel.fa

end=`date +%c`

# Mapping RNA-seq reads
gsnap -d /project/projectdirectory/mydirectory/dot_files/software/gmapdb/b73/ /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/*_1.fastq.gz /project/projectdirectory/mydirectory/rnaseq/maize/raw_data/*_2.fastq.gz
