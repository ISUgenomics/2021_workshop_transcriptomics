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



### Alignment rates
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="PAIR" {print $7}' ) |less
SRR1573504_1Gene.bam_alignment.stats    0.951596
SRR1573505_1Gene.bam_alignment.stats    0.933987
SRR1573506_1Gene.bam_alignment.stats    0.920971
SRR1573507_1Gene.bam_alignment.stats    0.925351
SRR1573508_1Gene.bam_alignment.stats    0.930736
SRR1573509_1Gene.bam_alignment.stats    0.948775
SRR1573510_1Gene.bam_alignment.stats    0.936508
SRR1573511_1Gene.bam_alignment.stats    0.917547
SRR1573512_1Gene.bam_alignment.stats    0.902367
SRR1573513_1Gene.bam_alignment.stats    0.945952
SRR1573514_1Gene.bam_alignment.stats    0.930853
SRR1573515_1Gene.bam_alignment.stats    0.937116
SRR1573516_1Gene.bam_alignment.stats    0.835353
SRR1573517_1Gene.bam_alignment.stats    0.896266
SRR1573518_1Gene.bam_alignment.stats    0.944063
SRR1573519_1Gene.bam_alignment.stats    0.942744
SRR1573520_1Gene.bam_alignment.stats    0.901736
SRR1573521_1Gene.bam_alignment.stats    0.933744
SRR1573522_1Gene.bam_alignment.stats    0.906409
SRR1573523_1Gene.bam_alignment.stats    0.943627
SRR1573524_1Gene.bam_alignment.stats    0.925833
SRR1573525_1Gene.bam_alignment.stats    0.92952
SRR1573526_1Gene.bam_alignment.stats    0.931899
SRR1573527_1Gene.bam_alignment.stats    0.933866
```

#### Assignment of unique reads to genes
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="PAIR" {print $2}' )  <(cat *Gene_counts_genes.txt.summary|awk '$1=="Status" {print $2}' ) <(cat *Gene_counts_genes.txt.summary | awk '$1=="Assigned" {print $2}' ) <(cat *Gene_counts_genes.txt
.summary | awk '$1=="Unassigned" {print $2}' ) |awk '{print $0,$4/$2}' |less

SRR1573504_1Gene.bam_alignment.stats    127954666       SRR1573504_1Gene_sorted.bam     113725482        0.888795
SRR1573505_1Gene.bam_alignment.stats    58423602        SRR1573505_1Gene_sorted.bam     50935698         0.871834
SRR1573506_1Gene.bam_alignment.stats    58765516        SRR1573506_1Gene_sorted.bam     50509286         0.859506
SRR1573507_1Gene.bam_alignment.stats    98029896        SRR1573507_1Gene_sorted.bam     84604160         0.863044
SRR1573508_1Gene.bam_alignment.stats    50082310        SRR1573508_1Gene_sorted.bam     43501035         0.868591
SRR1573509_1Gene.bam_alignment.stats    128128680       SRR1573509_1Gene_sorted.bam     112962496        0.881633
SRR1573510_1Gene.bam_alignment.stats    41381958        SRR1573510_1Gene_sorted.bam     36149951         0.873568
SRR1573511_1Gene.bam_alignment.stats    42928470        SRR1573511_1Gene_sorted.bam     36836120         0.858081
SRR1573512_1Gene.bam_alignment.stats    28356672        SRR1573512_1Gene_sorted.bam     23799585         0.839294
SRR1573513_1Gene.bam_alignment.stats    84059896        SRR1573513_1Gene_sorted.bam     74193550         0.882627
SRR1573514_1Gene.bam_alignment.stats    40178336        SRR1573514_1Gene_sorted.bam     34737709         0.864588
SRR1573515_1Gene.bam_alignment.stats    64372620        SRR1573515_1Gene_sorted.bam     55807294         0.866941
SRR1573516_1Gene.bam_alignment.stats    49097676        SRR1573516_1Gene_sorted.bam     38055327         0.775094
SRR1573517_1Gene.bam_alignment.stats    39684356        SRR1573517_1Gene_sorted.bam     32981690         0.831101
SRR1573518_1Gene.bam_alignment.stats    101234800       SRR1573518_1Gene_sorted.bam     88587226         0.875067
SRR1573519_1Gene.bam_alignment.stats    80364444        SRR1573519_1Gene_sorted.bam     70722826         0.880026
SRR1573520_1Gene.bam_alignment.stats    28985818        SRR1573520_1Gene_sorted.bam     20785283         0.717085
SRR1573521_1Gene.bam_alignment.stats    37739332        SRR1573521_1Gene_sorted.bam     32775917         0.868482
SRR1573522_1Gene.bam_alignment.stats    29619664        SRR1573522_1Gene_sorted.bam     24932950         0.84177
SRR1573523_1Gene.bam_alignment.stats    9449168 SRR1573523_1Gene_sorted.bam     8297034  0.87807
SRR1573524_1Gene.bam_alignment.stats    56449850        SRR1573524_1Gene_sorted.bam     47975143         0.849872
SRR1573525_1Gene.bam_alignment.stats    15903366        SRR1573525_1Gene_sorted.bam     13689548         0.860796
SRR1573526_1Gene.bam_alignment.stats    45958470        SRR1573526_1Gene_sorted.bam     39805049         0.866109
SRR1573527_1Gene.bam_alignment.stats    28339974        SRR1573527_1Gene_sorted.bam     24664579         0.870311
```

#### Assignment of unique and multiple mapping reads to genes
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="PAIR" {print $2}' )  <(cat *GeneMult_counts_genes.txt.summary|awk '$1=="Status" {print $2}' ) <(cat *GeneMult_counts_genes.txt.summary | awk '$1=="Assigned" {print $2}' ) <(cat *GeneMult_counts_genes.txt.summary | awk '$1=="Unassigned" {print $2}' ) |awk '{print $0,$4/$2}' |less

SRR1573504_1Gene.bam_alignment.stats    127954666       SRR1573504_1Gene_sorted.bam     119833441        0.93653
SRR1573505_1Gene.bam_alignment.stats    58423602        SRR1573505_1Gene_sorted.bam     53691682         0.919007
SRR1573506_1Gene.bam_alignment.stats    58765516        SRR1573506_1Gene_sorted.bam     53002346         0.901929
SRR1573507_1Gene.bam_alignment.stats    98029896        SRR1573507_1Gene_sorted.bam     88995926         0.907845
SRR1573508_1Gene.bam_alignment.stats    50082310        SRR1573508_1Gene_sorted.bam     45712896         0.912755
SRR1573509_1Gene.bam_alignment.stats    128128680       SRR1573509_1Gene_sorted.bam     119024010        0.928941
SRR1573510_1Gene.bam_alignment.stats    41381958        SRR1573510_1Gene_sorted.bam     38090018         0.92045
SRR1573511_1Gene.bam_alignment.stats    42928470        SRR1573511_1Gene_sorted.bam     38944593         0.907197
SRR1573512_1Gene.bam_alignment.stats    28356672        SRR1573512_1Gene_sorted.bam     25042392         0.883122
SRR1573513_1Gene.bam_alignment.stats    84059896        SRR1573513_1Gene_sorted.bam     78378885         0.932417
SRR1573514_1Gene.bam_alignment.stats    40178336        SRR1573514_1Gene_sorted.bam     36801046         0.915943
SRR1573515_1Gene.bam_alignment.stats    64372620        SRR1573515_1Gene_sorted.bam     59321058         0.921526
SRR1573516_1Gene.bam_alignment.stats    49097676        SRR1573516_1Gene_sorted.bam     40425815         0.823375
SRR1573517_1Gene.bam_alignment.stats    39684356        SRR1573517_1Gene_sorted.bam     35031851         0.882762
SRR1573518_1Gene.bam_alignment.stats    101234800       SRR1573518_1Gene_sorted.bam     94102329         0.929545
SRR1573519_1Gene.bam_alignment.stats    80364444        SRR1573519_1Gene_sorted.bam     74401025         0.925795
SRR1573520_1Gene.bam_alignment.stats    28985818        SRR1573520_1Gene_sorted.bam     31647525         1.09183
SRR1573521_1Gene.bam_alignment.stats    37739332        SRR1573521_1Gene_sorted.bam     35020472         0.927957
SRR1573522_1Gene.bam_alignment.stats    29619664        SRR1573522_1Gene_sorted.bam     26332280         0.889013
SRR1573523_1Gene.bam_alignment.stats    9449168 SRR1573523_1Gene_sorted.bam     8825328  0.933979
SRR1573524_1Gene.bam_alignment.stats    56449850        SRR1573524_1Gene_sorted.bam     52822743         0.935746
SRR1573525_1Gene.bam_alignment.stats    15903366        SRR1573525_1Gene_sorted.bam     14583785         0.917025
SRR1573526_1Gene.bam_alignment.stats    45958470        SRR1573526_1Gene_sorted.bam     42279641         0.919953
SRR1573527_1Gene.bam_alignment.stats    28339974        SRR1573527_1Gene_sorted.bam     26049721         0.919186

```
hmm.  one is above 100%, so is it counting multiple mapping reads for each feature they map to?



#### Assignment of unique reads to mRNAs
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="PAIR" {print $2}' )  <(cat *Mrna_counts_mrnas.txt.summary|awk '$1=="Status" {print $2}' ) <(cat *Mrna_counts_mrnas.txt.summary | awk '$1=="Assigned" {print $2}' ) <(cat *Mrna_counts_mrnas.txt.summary | awk '$1=="Unassigned" {print $2}' ) |awk '{print $0,$4/$2}' |less


SRR1573504_1Gene.bam_alignment.stats    127954666       SRR1573504_1Gene_sorted.bam     80538062         0.629427
SRR1573505_1Gene.bam_alignment.stats    58423602        SRR1573505_1Gene_sorted.bam     36377289         0.622647
SRR1573506_1Gene.bam_alignment.stats    58765516        SRR1573506_1Gene_sorted.bam     36009313         0.612763
SRR1573507_1Gene.bam_alignment.stats    98029896        SRR1573507_1Gene_sorted.bam     59919045         0.611232
SRR1573508_1Gene.bam_alignment.stats    50082310        SRR1573508_1Gene_sorted.bam     30949514         0.617973
SRR1573509_1Gene.bam_alignment.stats    128128680       SRR1573509_1Gene_sorted.bam     78162094         0.610028
SRR1573510_1Gene.bam_alignment.stats    41381958        SRR1573510_1Gene_sorted.bam     25894608         0.625746
SRR1573511_1Gene.bam_alignment.stats    42928470        SRR1573511_1Gene_sorted.bam     26123075         0.608526
SRR1573512_1Gene.bam_alignment.stats    28356672        SRR1573512_1Gene_sorted.bam     16830802         0.593539
SRR1573513_1Gene.bam_alignment.stats    84059896        SRR1573513_1Gene_sorted.bam     52417332         0.623571
SRR1573514_1Gene.bam_alignment.stats    40178336        SRR1573514_1Gene_sorted.bam     23433151         0.583229
SRR1573515_1Gene.bam_alignment.stats    64372620        SRR1573515_1Gene_sorted.bam     37426366         0.581402
SRR1573516_1Gene.bam_alignment.stats    49097676        SRR1573516_1Gene_sorted.bam     26361129         0.536912
SRR1573517_1Gene.bam_alignment.stats    39684356        SRR1573517_1Gene_sorted.bam     22401340         0.564488
SRR1573518_1Gene.bam_alignment.stats    101234800       SRR1573518_1Gene_sorted.bam     59871892         0.591416
SRR1573519_1Gene.bam_alignment.stats    80364444        SRR1573519_1Gene_sorted.bam     50870816         0.633002
SRR1573520_1Gene.bam_alignment.stats    28985818        SRR1573520_1Gene_sorted.bam     14462555         0.498953
SRR1573521_1Gene.bam_alignment.stats    37739332        SRR1573521_1Gene_sorted.bam     23463941         0.621737
SRR1573522_1Gene.bam_alignment.stats    29619664        SRR1573522_1Gene_sorted.bam     17619066         0.594844
SRR1573523_1Gene.bam_alignment.stats    9449168 SRR1573523_1Gene_sorted.bam     5908088  0.62525
SRR1573524_1Gene.bam_alignment.stats    56449850        SRR1573524_1Gene_sorted.bam     34241462         0.606582
SRR1573525_1Gene.bam_alignment.stats    15903366        SRR1573525_1Gene_sorted.bam     9657008  0.60723
SRR1573526_1Gene.bam_alignment.stats    45958470        SRR1573526_1Gene_sorted.bam     28319707         0.616202
SRR1573527_1Gene.bam_alignment.stats    28339974        SRR1573527_1Gene_sorted.bam     17697990         0.624489

```
