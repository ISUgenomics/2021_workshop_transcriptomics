# Transfer Maize and Bumblebee data to the Nova HPC cluster

1. Log onto Nova HPC

```sh
ssh kmcelroy@nova.its.iastate.edu
```

2. Make directories to work in

```sh
mkdir /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics
cd /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics
mkdir 00_genome
mkdir 01_data_maize
mkdir 02_data_bee
mkdir 03_qc_maize
mkdir 04_qc_bee
```

3. Download the Maize genome, transcriptome, and annotation (in GFF) from NCBI. Unzip each of these files. Kallisto and Salmon only need the transcriptome but I decided to download these other files that might be useful.

The Maize reference assembly used here is B73 Maize reference genome.

* NCBI genome page for Maize - [https://www.ncbi.nlm.nih.gov/assembly/GCF_902167145.1/](https://www.ncbi.nlm.nih.gov/genome/?term=txid4577[orgn])

```sh
# pwd = /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics/00_genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
gunzip GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
gunzip GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rna.fna.gz
gunzip GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rna.fna.gz
```

4. Same approach for collecting the Bumblebee reference assemblies.

* NCBI genome page for Bombus impatiens - [https://www.ncbi.nlm.nih.gov/genome/3415?genome_assembly_id=468201]

```sh
# pwd = /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics/00_genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_genomic.fna.gz
gunzip GCF_000188095.3_BIMP_2.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_rna.fna.gz
gunzip GCF_000188095.3_BIMP_2.2_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_genomic.gff.gz
gunzip GCF_000188095.3_BIMP_2.2_genomic.gff.gz
```

5. Collect the Maize RNA-seq data from NCBI SRA

Run this SLURM array script to download the SRA file with fastq-dump from the SRA-toolkit.

```sh
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 1 processor core(s) per node
#SBATCH --job-name="getSRA"
#SBATCH --mail-user=kmcelroy@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="getSRA" # job standard output file (%j replaced by job id)
#SBATCH --array=0-23 # run 24 jobs in an array

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load sra-toolkit/2.9.6-ub7kz5h

cd /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics/01_data

# tsv file downloaded from the workshop repo
cut -f 7 tsv | more | perl -pe 's/;/\n/g' | xargs -I xx echo "wget http://xx"
awk '{print $3}' tsv > sra

names=($(cat sra))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
fastq-dump --split-files --origfmt --gzip ${names[${SLURM_ARRAY_TASK_ID}]}
```

Run the above script:
```sh
cd /work/LAS/serb-lab/kmcelroy/2021_workshop_transcriptomics
sbatch getSRA.sbatch
```

6. Collect the Bumblebee RNA-seq data. I have the data but had to get it in a bit of a roundabout way ...
