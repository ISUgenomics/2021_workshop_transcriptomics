### Maize RNAseq:

Working directory: `/work/gif/Siva/Transcriptomic_Worshop`
* Directory Contents:

```
[csiva@novawide030 Transcriptomic_Worshop]$ ls -lth
.......

lrwxrwxrwx. 1 csiva its-hpc-nova-gif  44 Feb 16 12:13 Bumblebee_data -> /work/gif/TranscriptomicsWorkshop/bumbleBee/
lrwxrwxrwx. 1 csiva its-hpc-nova-gif  40 Feb 16 12:13 maize_data -> /work/gif/TranscriptomicsWorkshop/Maize/
```
* 1st Step: Quality check using fastqc

```
salloc -N1 -n36 -t 01:00:00
mkdir 01_QC
module load parallel
module load fastqc

parallel -j32 "fastqc -o 01_QC" ::: maize_data/*fastq

```

```
-rw-r--r--. 1 csiva its-hpc-nova-gif 42K Feb 16 13:05 fq.log
drwxr-sr-x. 2 csiva its-hpc-nova-gif  98 Feb 16 13:05 01_QC
```
fastqc makes a html and zipped file for each sample, in this case there are 48 files, so we end up with 96 files. Each html file has got several graphs. So we use multiqc to collate all the data.

#### Collating using multiqc
In Nova I used singularity to pull a docker image and then use that
```
cd 01_QC

singularity exec --bind $PWD /work/gif/Siva/multiqc_latest.sif multiqc .

multiqc : This is MultiQC v1.9
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching   : /work/gif/Siva/Transcriptomic_Worshop/01_QC
Searching 96 files..  [####################################]  100%          
[INFO   ]          fastqc : Found 48 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete


```
![multiqc.html](Figures/multiqc_report_maize_ligule.html)
