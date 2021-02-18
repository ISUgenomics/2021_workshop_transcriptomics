# Measure transcript abundance with Kallisto and Salmon for the Maize and bee data

First, make the indexes of each species transcriptome for Kallisto and Salmon

```sh
#load the modules for Kallisto and Salmon
module load kallisto/0.46.2-openmpi3-py3-bjls6kt
# need to load this gcc library for Salmon
module load gcc/7.3.0-xegsmw4
module load salmon/0.9.1-py2-vh2miiy
cd /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics/00_genome
kallisto index -i maize_transcripts.idx GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rna.fna
kallisto index -i bombus_transcripts.idx GCF_000188095.3_BIMP_2.2_rna.fna
# make an index with duplicate transcripts, to compare more directly to kallisto
salmon index --keepDuplicates -t GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rna.fna -i maize_index_dups
salmon index --keepDuplicates -t GCF_000188095.3_BIMP_2.2_rna.fna -i bombus_index_dups
```

## Maize data

1. Kallisto
Run this array job script to run kallisto quant on every pair of reads

```sh

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --job-name="kallisto_maize"
#SBATCH --mail-user=kmcelroy@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="kallisto_maize" # job standard output file (%j replaced by job id)
#SBATCH --array=0-23 # run 24 jobs in an array

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load kallisto/0.46.2-openmpi3-py3-bjls6kt

cd /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics

names=($(cat sra))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

kallisto quant -i 00_genome/maize_transcripts.idx \
        -o "05_kallisto_maize/${names[${SLURM_ARRAY_TASK_ID}]}" \
        -b 100 -t 16 --pseudobam \
        "03_qc_maize/fastp/${names[${SLURM_ARRAY_TASK_ID}]}_1.trim.fastq.gz" \
        "03_qc_maize/fastp/${names[${SLURM_ARRAY_TASK_ID}]}_2.trim.fastq.gz"
```

Run the above script:
```sh
/work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics
sbatch kallisto_maize.sbatch
```

2. Salmon
Run this array job script to run salmon quant on every pair of reads

```sh
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --job-name="salmon_maize"
#SBATCH --mail-user=kmcelroy@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="salmon_maize" # job standard output file (%j replaced by job id)
#SBATCH --array=0-23 # run 24 jobs in an array

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load gcc/7.3.0-xegsmw4
module load salmon/0.9.1-py2-vh2miiy

cd /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics
names=($(cat sra))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

salmon quant -i 00_genome/maize_index_dups -l A \
        -1 "03_qc_maize/fastp/${names[${SLURM_ARRAY_TASK_ID}]}_1.trim.fastq.gz" \
        -2 "03_qc_maize/fastp/${names[${SLURM_ARRAY_TASK_ID}]}_2.trim.fastq.gz" \
        -p 16 -o "06_salmon_maize/${names[${SLURM_ARRAY_TASK_ID}]}_quant"
```

Run the above script:
```sh
/work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics
sbatch salmon_maize.sbatch
```
