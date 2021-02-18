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

How similar are the Kallisto and Salmon abundance estimates? I generated one file with all of the counts and tpm values organized by SRA file, analysis program, and transcritp.

```sh
# make file with headers for columns
echo -e 'sample\ttranscript\tprogram\tcounts\ttpm' > tpm_comp.tab
# add all of the Kallisto output
for f in $(cat sra); do awk -v OFS='\t' 'FNR>1 {print var,$1,"Kallisto",$4,$5}' var=$f "05_kallisto_maize/$f/abundance.tsv" >> tpm_comp.tab; done
# add all of the Salmon output
for f in $(cat sra); do awk -v OFS='\t' 'FNR>1 {print var,$1,"Salmon",$5,$4}' var=$f "06_salmon_maize/${f}_quant/quant.sf" >> tpm_comp.tab; done

# grab an arbitrary (first) transcript from the above output file
# sample  transcript      program counts  tpm
grep NM_001111367.2 tpm_comp.tab | sort -k1
SRX699505       NM_001111367.2  Kallisto        9813    176.996
SRX699505       NM_001111367.2  Salmon  9816.000000     175.090656
SRX699506       NM_001111367.2  Kallisto        2946    114.576
SRX699506       NM_001111367.2  Salmon  2946.000000     113.187955
SRX699507       NM_001111367.2  Kallisto        2811    111.601
SRX699507       NM_001111367.2  Salmon  2811.000000     110.526127
SRX699508       NM_001111367.2  Kallisto        11086   269.265
SRX699508       NM_001111367.2  Salmon  11086.000000    266.657318
SRX699509       NM_001111367.2  Kallisto        3481.88 165.021
SRX699509       NM_001111367.2  Salmon  3481.879818     163.232900
SRX699510       NM_001111367.2  Kallisto        10926   214.428
SRX699510       NM_001111367.2  Salmon  10926.000000    212.311727
SRX699511       NM_001111367.2  Kallisto        4149    231.051
SRX699511       NM_001111367.2  Salmon  4149.000000     228.576818
SRX699512       NM_001111367.2  Kallisto        7021.96 392.848
SRX699512       NM_001111367.2  Salmon  7021.960761     388.493648
SRX699513       NM_001111367.2  Kallisto        2810    235.104
SRX699513       NM_001111367.2  Salmon  2810.000000     232.760154
SRX699514       NM_001111367.2  Kallisto        12327   366.291
SRX699514       NM_001111367.2  Salmon  12329.000000    362.852931
SRX699515       NM_001111367.2  Kallisto        3282    217.006
SRX699515       NM_001111367.2  Salmon  3282.000000     215.054077
SRX699516       NM_001111367.2  Kallisto        5828    239.781
SRX699516       NM_001111367.2  Salmon  5827.000000     237.676423
SRX699517       NM_001111367.2  Kallisto        3009    168.186
SRX699517       NM_001111367.2  Salmon  3009.000000     166.911977
SRX699518       NM_001111367.2  Kallisto        2740    188.504
SRX699518       NM_001111367.2  Salmon  2740.000000     187.021143
SRX699519       NM_001111367.2  Kallisto        8920    230.981
SRX699519       NM_001111367.2  Salmon  8922.000000     228.825612
SRX699520       NM_001111367.2  Kallisto        4265    127.3
SRX699520       NM_001111367.2  Salmon  4266.000000     126.096234
SRX699521       NM_001111367.2  Kallisto        1687    159.021
SRX699521       NM_001111367.2  Salmon  1687.000000     169.818033
SRX699522       NM_001111367.2  Kallisto        1277    80.9199
SRX699522       NM_001111367.2  Salmon  1277.000000     80.232599
SRX699523       NM_001111367.2  Kallisto        1935    166.133
SRX699523       NM_001111367.2  Salmon  1935.000000     164.446554
SRX699524       NM_001111367.2  Kallisto        251     63.4539
SRX699524       NM_001111367.2  Salmon  251.000000      62.885420
SRX699525       NM_001111367.2  Kallisto        3600    156.271
SRX699525       NM_001111367.2  Salmon  3600.000000     155.837192
SRX699526       NM_001111367.2  Kallisto        2180    338.562
SRX699526       NM_001111367.2  Salmon  2180.000000     335.637095
SRX699527       NM_001111367.2  Kallisto        6054    318.93
SRX699527       NM_001111367.2  Salmon  6054.000000     316.118090
SRX699528       NM_001111367.2  Kallisto        5885    504.337
SRX699528       NM_001111367.2  Salmon  5885.000000     499.388370
```
