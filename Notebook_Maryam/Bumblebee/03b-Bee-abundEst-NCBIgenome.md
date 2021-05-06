# Abundance estimation

I need to rerun with NCBI genome to be able to compare with RIck and Jennifer data :


### converting the sam files to bam files
First we need to convert the alignment sam files to bam files:


```bash
module load samtools
for n in ../02Bee-STAR/*.sam ; do   m=`echo $n|  cut -d'/' -f 3 | cut -d'-' -f 3,4  | cut -d'.' -f 1`; samtools view --threads 12 -bS -o ./bam/$m.bam $n; done

```
### Count estimate

Now that I have the bam files I can proceed with `featureCounts`.

```
module load subread
module load parallel

parallel -j 12 "featureCounts -T 4 -s 0 -p -t gene -g ID -a /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/Bombus_impatiens_BIMP_2.2_RefSeq_proteincoding.gff3 -o ./counts/{/.}.counts.txt {}" ::: bam/*.bam
```

options :
```
-T    Number of the threads
-s    perform strand-specific read counts. 0 for single end
-t  feature type in GTF
-g   attribute type in GTF
```

### combining count files

extracting the count column from each count file:

```
for n in *.txt ; do  m=`echo $n|  cut -d'-' -f 1` ; awk 'BEGIN {OFS="\t"} {print $7}' $n > t-$m.txt ; done
```

combine all new files with the ID to generate a single table of counts :

```
paste <(awk 'BEGIN {OFS="\t"} {print $1}' E12-F10Aligned.counts.txt) t-*.txt | grep -v '^\#' > At_count.txt
```
