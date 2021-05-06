# Abundance estimation

* Nova: /work/gif/Maryam/projects/Transcriptomics/03Maize-abunEst
* 22Feb , 2021

### converting the sam files to bam files
First we need to convert the alignment sam files to bam files:


```bash
salloc -N1 -n12 -p DEBUG -t 1:00:00

cd /work/gif/Maryam/projects/Transcriptomics/03Maize-abunEst

module load samtools
for n in ../02Maize-STAR/StarOut*.sam ; do   m=`echo $n|  cut -d'/' -f 3 | cut -d'-' -f 3,4  | cut -d'.' -f 1`; samtools view --threads 12 -bS -o ./bam/$m.bam $n; done

```

### Count estimate

Now that I have the bam files I can proceed with `featureCounts`.

```
module load subread
module load parallel

parallel -j 12 "featureCounts -T 4 -s 1 -p -t gene -g ID -a /work/gif/Maryam/projects/Transcriptomics/00-rawdata/maize/Bombus_impatiens_BIMP_2.2_RefSeq_proteincoding.gff3 -o ./counts/{/.}.counts.txt {}" ::: bam/*.bam
```
