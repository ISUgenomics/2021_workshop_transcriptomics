# Align maize leaf sheath data to genome


### prepare for alignment
```
#/work/gif/TranscriptomicsWorkshop/remkv6/01_Maize

# gather genome and fastq
ln -s /work/gif/TranscriptomicsWorkshop/Maize/01_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna Maize.fa
ln -s /work/gif/TranscriptomicsWorkshop/Maize/01_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff
ml hisat2; hisat2-build Maize.fa Maize

for f in /work/gif/TranscriptomicsWorkshop/Maize/*fastq; do ln -s $f;done


```

### Alignment script
```
#runFeatureCountsPairedUnstranded.sh
################################################################################################################################
#note premake the database, hisat2-build genome.fa genome
#Script expects the database to be named the same as the genome file without the last extension (Genome), not (Genome.fa)


#Example run

#sh runFeatureCountsPairedUnstranded.sh {read1.fastq} {read2.fastq} {Directory of Database, not actual database without / at the end} {genomeName} {Genic GFF}

PROC=36
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"
GFF="$5"


module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p ${PROC}  -x ${GENOME%.*} -1 ${R1_FQ} -2 ${R2_FQ}  -S ${R1_FQ%.*}Gene.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}Gene.bam ${R1_FQ%.*}Gene.sam
mkdir ${R1_FQ%.*}Gene_temp
samtools sort  -o ${R1_FQ%.*}Gene_sorted.bam -T ${R1_FQ%.*}Gene_temp --threads ${PROC} ${R1_FQ%.*}Gene.bam


module load picard/2.17.0-ft5qztz
java -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/picard-2.17.0-ft5qztzntoymuxiqt3b6yi6uqcmgzmds/bin/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${R1_FQ%.*}Gene_sorted.bam OUTPUT=${R1_FQ%.*}Gene.bam_alignment.stats

module load subread
featureCounts -T ${PROC}   -t gene -g ID -a ${GFF} -o ${R1_FQ%.*}Gene_counts_genes.txt ${R1_FQ%.*}Gene_sorted.bam
featureCounts -T ${PROC}   -t mRNA  -g ID -a ${GFF} -o ${R1_FQ%.*}Mrna_counts_mrnas.txt ${R1_FQ%.*}Gene_sorted.bam
featureCounts -T ${PROC} -M   -t gene -g ID -a ${GFF} -o ${R1_FQ%.*}GeneMult_counts_genes.txt ${R1_FQ%.*}Gene_sorted.bam

################################################################################################################################
```

### Align
```
paste <(ls -1 *_1.fastq) <(ls -1 *_2.fastq) |sed 's/\t/ /g' |while read line ; do echo "sh runFeatureCountsPairedUnstranded.sh "$line" /work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align Maize.fa GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff";done >Maizealign.sh

```
