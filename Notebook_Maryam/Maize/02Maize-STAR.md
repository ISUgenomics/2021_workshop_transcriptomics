# Alignment with Star Maize dataset
* Feb 2021
* Nova : /work/gif/Maryam/projects/Transcriptomics/02Maize-STAR

### Ref data set:
* ref genome : `GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna`


#### First step: Generating genome index
Building the index:

```bash
salloc -N1 -n8 -p DEBUG -t 1:00:00
/work/gif/Maryam/projects/Transcriptomics/02Maize-STAR
module load star

mkdir MaizeDB

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./MaizeDB  --sjdbGTFfile /work/gif/Maryam/projects/Transcriptomics/00-rawdata/maize/01_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff  --genomeFastaFiles /work/gif/Maryam/projects/Transcriptomics/00-rawdata//maize/01_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna --sjdbOverhang  99

```
* Note: For details on the parameters visit `02-Bee-STAR.md`.

#### Second Step: Alignment

```bash
#!/bin/bash

read1="$1"

filename=`echo $read1 | cut -d'/' -f 4|  cut -d'_' -f 1`
jobname=$filename
filename+=".sub"



echo $filename
echo "file name $filename"

echo "#!/bin/bash" > $filename
echo "#SBATCH -N 1" >> $filename
echo "#SBATCH --ntasks-per-node=12" >> $filename
echo "#SBATCH -t 24:00:00" >>  $filename
echo "#SBATCH -J star-$jobname " >> $filename
echo "#SBATCH -o ./log/star-$jobname.o%j  " >> $filename
echo "#SBATCH -e ./log/star-$jobname.e%j " >> $filename
echo "#SBATCH --mail-user=msayadi@iastate.edu  " >> $filename
echo "#SBATCH --mail-type=begin " >> $filename
echo "#SBATCH --mail-type=end " >> $filename

echo >> $filename



echo "cd \$SLURM_SUBMIT_DIR" >> $filename
echo "ulimit -s unlimited" >> $filename
echo "scontrol show job \$SLURM_JOB_ID" >> $filename

read2=$(echo $read1 | sed 's/_1/_2/g');

echo >> $filename
echo "module load star"  >> $filename
echo >> $filename

echo "star --runThreadN 12 --readFilesIn $read1 $read2  --genomeDir /work/gif/Maryam/projects/Transcriptomics/02Maize-STAR/MaizeDB  --outFileNamePrefix Star" >> $filename
echo "scontrol show job \$SLURM_JOB_ID" >> $filename
```

All the nodes were busy on Nova so I decided to run them on an interactive node:

```bash
for n in ../00-rawdata/maize/*.fastq; do   read1=$n; filename=`echo $read1 | cut -d'/' -f 4|  cut -d'_' -f 1`; read2=`echo $read1 | sed 's/_1/_2/g'`; STAR --runThreadN 24  --readFilesIn $read1 $read2 --genomeDir /work/gif/Maryam/projects/Transcriptomics/02Maize-STAR/MaizeDB   --outFileNamePrefix  Star-$filename; done

```
