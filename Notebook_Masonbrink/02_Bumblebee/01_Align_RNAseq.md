# Gather and align bumblebee data
```
#/work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/01_Align

#softlink and make genome database
ln -s /work/gif/TranscriptomicsWorkshop/bumbleBee/01_Genome/GCF_000188095.3_BIMP_2.2_genomic.fna Bumblebee.fasta
ml hisat2; hisat2-build Bumblebee.fasta Bumblebee


#Softlink fastq files
for f in /work/gif/TranscriptomicsWorkshop/bumbleBee/*gz ; do ln -s $f;done


cp /work/gif/TranscriptomicsWorkshop/bumbleBee/01_Genome/GCF_000188095.3_BIMP_2.2_genomic.gff.gz .
gunzip GCF_000188095.3_BIMP_2.2_genomic.gff.gz
```


### Alignment script
```
#!/bin/bash
#runFeatureCountsSingleEnd.sh
################################################################################################################################
#note premake the database, hisat2-build genome.fa genome
#Script expects the database to be named the same as the genome file without the last extension (Genome), not (Genome.fa)


#Example run

#sh runFeatureCountsSingleEnd.sh {readsFile} {Directory of Database, not actual database without / at the end} {genomeName} {Genic GFF}


PROC=36
READS_FQ="$1"
DBDIR="$2"
GENOME="$3"
GFF="$4"


module load hisat2
hisat2 -p ${PROC}  -x ${GENOME%.*} -U ${READS_FQ}  -S ${READS_FQ%.*}Gene.sam

module load samtools
samtools view --threads ${PROC} -b -o ${READS_FQ%.*}Gene.bam ${READS_FQ%.*}Gene.sam
mkdir ${READS_FQ%.*}Gene_temp
samtools sort  -o ${READS_FQ%.*}Gene_sorted.bam -T ${READS_FQ%.*}Gene_temp --threads ${PROC} ${READS_FQ%.*}Gene.bam


module load picard/2.17.0-ft5qztz
java -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/picard-2.17.0-ft5qztzntoymuxiqt3b6yi6uqcmgzmds/bin/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${READS_FQ%.*}Gene_sorted.bam OUTPUT=${READS_FQ%.*}Gene.bam_alignment.stats

module load subread
featureCounts -T ${PROC}  -t gene  -g ID -a ${GFF} -o ${READS_FQ%.*}Gene_counts_genes.txt ${READS_FQ%.*}Gene_sorted.bam
featureCounts -T ${PROC}  -t mRNA  -g ID -a ${GFF} -o ${READS_FQ%.*}mRNA_counts_genes.txt ${READS_FQ%.*}Gene_sorted.bam
featureCounts -T ${PROC} -M   -t gene  -g ID -a ${GFF} -o ${READS_FQ%.*}GeneMult_counts_genes.txt ${READS_FQ%.*}Gene_sorted.bam

################################################################################################################################
```


### Align
```
for f in *gz; do echo "sh runFeatureCountsSingleEnd.sh "$f" /work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/01_Align Bumblebee.fasta GCF_000188095.3_BIMP_2.2_genomic.gff";done >align.sh
```
