2/9/21

bumblebee mapping

made hisat2 database building script
nano hisatbuild.sub

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -J HI_build
#SBATCH -o HI_build.o%j
#SBATCH -e HI_build.e%j
#SBATCH --mail-user=awalton@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
set -o xtrace
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited

module load hisat2

GENOME="/work/LAS/amytoth-lab/bumblebee/genome/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa"
hisat2-build $GENOME ${GENOME%.*}
scontrol show job $SLURM_JOB_ID

```
HiSat2 mapping script
file name  = hisatmap

```
#!/bin/bash
set -o xtrace
# set the reference index:
GENOME="work/LAS/amytoth-lab/bumblebee/genome/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa"
# make an output directory to store the output aligned files
mkdir -p hisatOut
# set that as the output directory
ODIR="hisatOut"


p=8 # use 8 threads
R1_FQ="$1" # first argument
R2_FQ="$2" # second argument

# purge and load relevant modules.
module purge

module load hisat2
module load samtools

OUTPUT=$(basename ${R1_FQ} |cut -f 1,2 -d "_");

hisat2 \
  -p ${p} \
  -x ${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ} \
  -S $ODIR\/${OUTPUT}.sam &> ${OUTPUT}.log
samtools view --threads 8 -bS -o $ODIR\/${OUTPUT}.bam $ODIR\/${OUTPUT}.sam

rm $ODIR\/${OUTPUT}.sam
```
