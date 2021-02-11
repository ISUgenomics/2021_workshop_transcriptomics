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

```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --job-name="hisatloop"
#SBATCH --mail-user=rfortune@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL



module purge
module load hisat2/2.2.0-5kvb7f2
module load samtools/1.10-py3-xuj7ylj
cd /work/LAS/amytoth-lab/awalton/transciptomics_class/01_bumble_HiSat2/bumblebee

for i in $(ls reads/*_L002_R1_001.fastq | uniq)
do
	hisat2   -p 16 -x genome/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic -U ${i} -S "SAMs${i}.sam";
done >& hisatloop.log
```
