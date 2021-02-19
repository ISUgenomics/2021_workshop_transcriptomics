02_feature featureCounts

##2-18-21
##working directory: /work/LAS/amytoth-lab/awalton/transciptomics_class/01_bumble_HiSat2/bumblebee/BAMs

#featurecounts

```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=2:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --job-name="featcount"
#SBATCH --mail-user=awalton@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load subread/1.6.0-ak6vxhs
cd /work/LAS/amytoth-lab/awalton/transciptomics_class/01_bumble_HiSat2/bumblebee/BAMs
for i in $(ls *.bam)
do
        featureCounts -T 16 -t gene -g ID -a /work/LAS/amytoth-lab/awalton/transciptomics_class/01_bumble_HiSat2/bumblebee/BAMs/Bombus_impatiens_BIMP_2.2_RefSeq_proteincoding.gff3 -o $i.counts.txt $i;
done
```

# Ryan's code for making a tidy counts filtered_1703
```
ut -f 1 1-A01-A1_S7_L002_R1_001.fastq.sam.bam.counts.txt | grep -v '#' > col1.txt
for i in *.counts.txt; do cat $i | grep -v '#' | cut -f 7 > $i.cut.txt ; done
paste col1.txt *cut* > count_table.txt
```
## and you don't need to sbatch these
