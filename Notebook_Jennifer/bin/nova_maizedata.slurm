#! /usr/bin/env bash
#SBATCH --nodes=1
# -p service
#SBATCH --ntasks-per-node=2
#SBATCH --time=24:00:00
#SBATCH --job-name=Maize
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
# --mail-user=someone@email.com
# --mail-type=begin
# --mail-type=end
# --account=projectname

set -e
set -u

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Main Program
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573506/SRR1573506_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573506/SRR1573506_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573507/SRR1573507_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573507/SRR1573507_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573508/SRR1573508_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573508/SRR1573508_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573509/SRR1573509_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573509/SRR1573509_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573510/SRR1573510_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573510/SRR1573510_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573511/SRR1573511_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573511/SRR1573511_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573512/SRR1573512_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573512/SRR1573512_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573513/SRR1573513_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573513/SRR1573513_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573514/SRR1573514_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573514/SRR1573514_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573515/SRR1573515_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573515/SRR1573515_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573516/SRR1573516_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573516/SRR1573516_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573517/SRR1573517_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573517/SRR1573517_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573518/SRR1573518_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573518/SRR1573518_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573519/SRR1573519_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573519/SRR1573519_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573520/SRR1573520_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573520/SRR1573520_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573521/SRR1573521_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573521/SRR1573521_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573522/SRR1573522_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573522/SRR1573522_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573523/SRR1573523_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573523/SRR1573523_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573524/SRR1573524_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573524/SRR1573524_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573525/SRR1573525_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573525/SRR1573525_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573526/SRR1573526_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573526/SRR1573526_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573527/SRR1573527_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573527/SRR1573527_2.fastq.gz
