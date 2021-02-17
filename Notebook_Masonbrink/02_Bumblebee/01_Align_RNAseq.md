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

#### Alignment stats
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="UNPAIRED" {print $7}' ) |less

1-A01-A1_S7_L002_R1_001.fastqGene.bam_alignment.stats   0.716843
1-A02-A2_S8_L002_R1_001.fastqGene.bam_alignment.stats   0.747593
1-A03-A3_S9_L002_R1_001.fastqGene.bam_alignment.stats   0.766884
1-A04-A4_S10_L002_R1_001.fastqGene.bam_alignment.stats  0.671936
1-A05-A5_S11_L002_R1_001.fastqGene.bam_alignment.stats  0.733566
1-A06-A6_S12_L002_R1_001.fastqGene.bam_alignment.stats  0.722409
1-A07-A7_S13_L002_R1_001.fastqGene.bam_alignment.stats  0.742856
1-A08-A8_S14_L002_R1_001.fastqGene.bam_alignment.stats  0.706004
1-A09-A9_S15_L002_R1_001.fastqGene.bam_alignment.stats  0.76727
1-A10-A10_S16_L002_R1_001.fastqGene.bam_alignment.stats 0.745161
1-A11-A11_S17_L002_R1_001.fastqGene.bam_alignment.stats 0.722867
1-A12-A12_S18_L002_R1_001.fastqGene.bam_alignment.stats 0.695275
1-B01-A13_S19_L002_R1_001.fastqGene.bam_alignment.stats 0.698609
1-B02-A14_S20_L002_R1_001.fastqGene.bam_alignment.stats 0.671933
1-B03-A16_S21_L002_R1_001-024.fastqGene.bam_alignment.stats     0.651383
1-B04-A17_S22_L002_R1_001.fastqGene.bam_alignment.stats 0.664987
1-B05-A18_S23_L002_R1_001.fastqGene.bam_alignment.stats 0.669463
1-B06-A19_S24_L002_R1_001.fastqGene.bam_alignment.stats 0.713439
1-B07-B2_S25_L002_R1_001.fastqGene.bam_alignment.stats  0.714408
1-B08-B3_S26_L002_R1_001.fastqGene.bam_alignment.stats  0.688258
1-B09-B4_S27_L002_R1_001.fastqGene.bam_alignment.stats  0.739815
1-B10-B5_S28_L002_R1_001.fastqGene.bam_alignment.stats  0.691651
1-B11-B6_S29_L002_R1_001.fastqGene.bam_alignment.stats  0.498513
1-B12-B7_S30_L002_R1_001.fastqGene.bam_alignment.stats  0.746438
1-C01-B8_S31_L002_R1_001.fastqGene.bam_alignment.stats  0.787653
1-C02-B10_S32_L002_R1_001.fastqGene.bam_alignment.stats 0.710902
1-C03-C1_S33_L002_R1_001.fastqGene.bam_alignment.stats  0.644792
1-C04-C2_S34_L002_R1_001.fastqGene.bam_alignment.stats  0.700882
1-C05-C3_S35_L002_R1_001.fastqGene.bam_alignment.stats  0.78088
1-C06-C4_S36_L002_R1_001.fastqGene.bam_alignment.stats  0.782866
1-C07-C6_S37_L002_R1_001.fastqGene.bam_alignment.stats  0.650172
1-C08-C7_S38_L002_R1_001.fastqGene.bam_alignment.stats  0.700423
1-C09-C8_S39_L002_R1_001.fastqGene.bam_alignment.stats  0.709526
1-C10-C9_S40_L002_R1_001.fastqGene.bam_alignment.stats  0.68879
1-C11-D1_S41_L002_R1_001.fastqGene.bam_alignment.stats  0.743857
1-C12-D2_S42_L002_R1_001.fastqGene.bam_alignment.stats  0.718843
1-D01-D5_S43_L002_R1_001.fastqGene.bam_alignment.stats  0.652962
1-D02-D6_S44_L002_R1_001.fastqGene.bam_alignment.stats  0.706298
1-D03-D7_S45_L002_R1_001.fastqGene.bam_alignment.stats  0.582135
1-D04-D8_S46_L002_R1_001.fastqGene.bam_alignment.stats  0.583949
1-D05-D9_S47_L002_R1_001.fastqGene.bam_alignment.stats  0.674012
1-D06-D10_S48_L002_R1_001.fastqGene.bam_alignment.stats 0.614037
1-D07-D11_S49_L002_R1_001.fastqGene.bam_alignment.stats 0.70471
1-D08-E1_S50_L002_R1_001.fastqGene.bam_alignment.stats  0.607073
1-D09-E2_S51_L002_R1_001.fastqGene.bam_alignment.stats  0.692269
1-D10-E5_S52_L002_R1_001.fastqGene.bam_alignment.stats  0.606752
1-D11-E7_S53_L002_R1_001.fastqGene.bam_alignment.stats  0.744952
1-D12-E8_S54_L002_R1_001.fastqGene.bam_alignment.stats  0.772314
1-E01-E9_S55_L002_R1_001.fastqGene.bam_alignment.stats  0.754534
1-E02-E10_S56_L002_R1_001.fastqGene.bam_alignment.stats 0.658656
1-E03-F1_S57_L002_R1_001.fastqGene.bam_alignment.stats  0.754044
1-E04-F2_S58_L002_R1_001.fastqGene.bam_alignment.stats  0.661595
1-E05-F3_S59_L002_R1_001.fastqGene.bam_alignment.stats  0.734738
1-E06-F4_S60_L002_R1_001.fastqGene.bam_alignment.stats  0.640634
1-E07-F5_S61_L002_R1_001.fastqGene.bam_alignment.stats  0.503371
1-E08-F6_S62_L002_R1_001.fastqGene.bam_alignment.stats  0.704253
1-E09-F7_S63_L002_R1_001.fastqGene.bam_alignment.stats  0.709601
1-E10-F8_S64_L002_R1_001.fastqGene.bam_alignment.stats  0.688254
1-E11-F9_S65_L002_R1_001.fastqGene.bam_alignment.stats  0.711303
1-E12-F10_S66_L002_R1_001.fastqGene.bam_alignment.stats 0.660544
```
Hmm.  50-78% alignment rate?  Bacterial contamnation?

### How many of the total unique mapped reads were assigned to a gene feature?
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="UNPAIRED" {print $2}' )  <(cat *Gene_counts_genes.txt.summary|awk '$1=="Status" {print $2}' ) <(cat *Gene_counts_genes.txt.summary | awk '$1=="Assigned" {print $2}' ) <(cat *Gene_counts_genes.txt.summary | awk '$1=="Unassigned" {print $2}' ) |awk '{print $0,$4/$2}' |less

1-A01-A1_S7_L002_R1_001.fastqGene.bam_alignment.stats   6103950 1-A01-A1_S7_L002_R1_001.fastqGene_sorted.bam    3112246  0.509874
1-A02-A2_S8_L002_R1_001.fastqGene.bam_alignment.stats   3800050 1-A02-A2_S8_L002_R1_001.fastqGene_sorted.bam    2052610  0.540153
1-A03-A3_S9_L002_R1_001.fastqGene.bam_alignment.stats   3724472 1-A03-A3_S9_L002_R1_001.fastqGene_sorted.bam    2014663  0.540926
1-A04-A4_S10_L002_R1_001.fastqGene.bam_alignment.stats  4604921 1-A04-A4_S10_L002_R1_001.fastqGene_sorted.bam   2238346  0.486077
1-A05-A5_S11_L002_R1_001.fastqGene.bam_alignment.stats  3510078 1-A05-A5_S11_L002_R1_001.fastqGene_sorted.bam   1836458  0.523196
1-A06-A6_S12_L002_R1_001.fastqGene.bam_alignment.stats  4919113 1-A06-A6_S12_L002_R1_001.fastqGene_sorted.bam   2562364  0.5209
1-A07-A7_S13_L002_R1_001.fastqGene.bam_alignment.stats  5920091 1-A07-A7_S13_L002_R1_001.fastqGene_sorted.bam   3140080  0.530411
1-A08-A8_S14_L002_R1_001.fastqGene.bam_alignment.stats  6352631 1-A08-A8_S14_L002_R1_001.fastqGene_sorted.bam   3209302  0.505193
1-A09-A9_S15_L002_R1_001.fastqGene.bam_alignment.stats  6078054 1-A09-A9_S15_L002_R1_001.fastqGene_sorted.bam   3357422  0.552384
1-A10-A10_S16_L002_R1_001.fastqGene.bam_alignment.stats 5341438 1-A10-A10_S16_L002_R1_001.fastqGene_sorted.bam  2878252  0.538853
1-A11-A11_S17_L002_R1_001.fastqGene.bam_alignment.stats 5527935 1-A11-A11_S17_L002_R1_001.fastqGene_sorted.bam  2864149  0.518123
1-A12-A12_S18_L002_R1_001.fastqGene.bam_alignment.stats 5515285 1-A12-A12_S18_L002_R1_001.fastqGene_sorted.bam  2717232  0.492673
1-B01-A13_S19_L002_R1_001.fastqGene.bam_alignment.stats 6509294 1-B01-A13_S19_L002_R1_001.fastqGene_sorted.bam  3276027  0.503285
1-B02-A14_S20_L002_R1_001.fastqGene.bam_alignment.stats 3897336 1-B02-A14_S20_L002_R1_001.fastqGene_sorted.bam  1890095  0.484971
1-B03-A16_S21_L002_R1_001-024.fastqGene.bam_alignment.stats     23324186        1-B03-A16_S21_L002_R1_001-024.fastqGene_sorted.bam      11017661         0.472371
1-B04-A17_S22_L002_R1_001.fastqGene.bam_alignment.stats 3598206 1-B04-A17_S22_L002_R1_001.fastqGene_sorted.bam  1733632  0.481805
1-B05-A18_S23_L002_R1_001.fastqGene.bam_alignment.stats 2951360 1-B05-A18_S23_L002_R1_001.fastqGene_sorted.bam  1425829  0.483109
1-B06-A19_S24_L002_R1_001.fastqGene.bam_alignment.stats 5570408 1-B06-A19_S24_L002_R1_001.fastqGene_sorted.bam  2876025  0.516304
1-B07-B2_S25_L002_R1_001.fastqGene.bam_alignment.stats  6208488 1-B07-B2_S25_L002_R1_001.fastqGene_sorted.bam   3192068  0.514146
1-B08-B3_S26_L002_R1_001.fastqGene.bam_alignment.stats  5910328 1-B08-B3_S26_L002_R1_001.fastqGene_sorted.bam   2894577  0.489749
1-B09-B4_S27_L002_R1_001.fastqGene.bam_alignment.stats  4150147 1-B09-B4_S27_L002_R1_001.fastqGene_sorted.bam   2194907  0.528875
1-B10-B5_S28_L002_R1_001.fastqGene.bam_alignment.stats  5304293 1-B10-B5_S28_L002_R1_001.fastqGene_sorted.bam   2593523  0.488948
1-B11-B6_S29_L002_R1_001.fastqGene.bam_alignment.stats  2544152 1-B11-B6_S29_L002_R1_001.fastqGene_sorted.bam   908769   0.357199
1-B12-B7_S30_L002_R1_001.fastqGene.bam_alignment.stats  7076685 1-B12-B7_S30_L002_R1_001.fastqGene_sorted.bam   3671100  0.51876
1-C01-B8_S31_L002_R1_001.fastqGene.bam_alignment.stats  5274379 1-C01-B8_S31_L002_R1_001.fastqGene_sorted.bam   2905700  0.550908
1-C02-B10_S32_L002_R1_001.fastqGene.bam_alignment.stats 4254713 1-C02-B10_S32_L002_R1_001.fastqGene_sorted.bam  2126398  0.499775
1-C03-C1_S33_L002_R1_001.fastqGene.bam_alignment.stats  4210123 1-C03-C1_S33_L002_R1_001.fastqGene_sorted.bam   1912824  0.454339
1-C04-C2_S34_L002_R1_001.fastqGene.bam_alignment.stats  2958147 1-C04-C2_S34_L002_R1_001.fastqGene_sorted.bam   1478037  0.49965
1-C05-C3_S35_L002_R1_001.fastqGene.bam_alignment.stats  4199024 1-C05-C3_S35_L002_R1_001.fastqGene_sorted.bam   2285592  0.544315
1-C06-C4_S36_L002_R1_001.fastqGene.bam_alignment.stats  4783587 1-C06-C4_S36_L002_R1_001.fastqGene_sorted.bam   2631494  0.550109
1-C07-C6_S37_L002_R1_001.fastqGene.bam_alignment.stats  4191393 1-C07-C6_S37_L002_R1_001.fastqGene_sorted.bam   1947280  0.46459
1-C08-C7_S38_L002_R1_001.fastqGene.bam_alignment.stats  4940838 1-C08-C7_S38_L002_R1_001.fastqGene_sorted.bam   2405608  0.486883
1-C09-C8_S39_L002_R1_001.fastqGene.bam_alignment.stats  4717779 1-C09-C8_S39_L002_R1_001.fastqGene_sorted.bam   2377656  0.503978
1-C10-C9_S40_L002_R1_001.fastqGene.bam_alignment.stats  5864332 1-C10-C9_S40_L002_R1_001.fastqGene_sorted.bam   2868607  0.489162
1-C11-D1_S41_L002_R1_001.fastqGene.bam_alignment.stats  5446043 1-C11-D1_S41_L002_R1_001.fastqGene_sorted.bam   2875324  0.527966
1-C12-D2_S42_L002_R1_001.fastqGene.bam_alignment.stats  4443892 1-C12-D2_S42_L002_R1_001.fastqGene_sorted.bam   2282615  0.513652
1-D01-D5_S43_L002_R1_001.fastqGene.bam_alignment.stats  5688965 1-D01-D5_S43_L002_R1_001.fastqGene_sorted.bam   2634357  0.463064
1-D02-D6_S44_L002_R1_001.fastqGene.bam_alignment.stats  3583246 1-D02-D6_S44_L002_R1_001.fastqGene_sorted.bam   1818758  0.507573
1-D03-D7_S45_L002_R1_001.fastqGene.bam_alignment.stats  4803348 1-D03-D7_S45_L002_R1_001.fastqGene_sorted.bam   1990017  0.414298
1-D04-D8_S46_L002_R1_001.fastqGene.bam_alignment.stats  3265963 1-D04-D8_S46_L002_R1_001.fastqGene_sorted.bam   1344631  0.41171
1-D05-D9_S47_L002_R1_001.fastqGene.bam_alignment.stats  4167619 1-D05-D9_S47_L002_R1_001.fastqGene_sorted.bam   1992712  0.478142
1-D06-D10_S48_L002_R1_001.fastqGene.bam_alignment.stats 4350687 1-D06-D10_S48_L002_R1_001.fastqGene_sorted.bam  1909747  0.438953
1-D07-D11_S49_L002_R1_001.fastqGene.bam_alignment.stats 8144167 1-D07-D11_S49_L002_R1_001.fastqGene_sorted.bam  4148343  0.509364
1-D08-E1_S50_L002_R1_001.fastqGene.bam_alignment.stats  4785525 1-D08-E1_S50_L002_R1_001.fastqGene_sorted.bam   2073003  0.433182
1-D09-E2_S51_L002_R1_001.fastqGene.bam_alignment.stats  4628678 1-D09-E2_S51_L002_R1_001.fastqGene_sorted.bam   2303411  0.497639
1-D10-E5_S52_L002_R1_001.fastqGene.bam_alignment.stats  4738622 1-D10-E5_S52_L002_R1_001.fastqGene_sorted.bam   2069042  0.436634
1-D11-E7_S53_L002_R1_001.fastqGene.bam_alignment.stats  2707274 1-D11-E7_S53_L002_R1_001.fastqGene_sorted.bam   1429658  0.52808
1-D12-E8_S54_L002_R1_001.fastqGene.bam_alignment.stats  5761711 1-D12-E8_S54_L002_R1_001.fastqGene_sorted.bam   3144512  0.54576
1-E01-E9_S55_L002_R1_001.fastqGene.bam_alignment.stats  3718478 1-E01-E9_S55_L002_R1_001.fastqGene_sorted.bam   1991105  0.535462
1-E02-E10_S56_L002_R1_001.fastqGene.bam_alignment.stats 5871181 1-E02-E10_S56_L002_R1_001.fastqGene_sorted.bam  2756625  0.469518
1-E03-F1_S57_L002_R1_001.fastqGene.bam_alignment.stats  5837908 1-E03-F1_S57_L002_R1_001.fastqGene_sorted.bam   3141162  0.538063
1-E04-F2_S58_L002_R1_001.fastqGene.bam_alignment.stats  10028269        1-E04-F2_S58_L002_R1_001.fastqGene_sorted.bam   4752280  0.473888
1-E05-F3_S59_L002_R1_001.fastqGene.bam_alignment.stats  2743201 1-E05-F3_S59_L002_R1_001.fastqGene_sorted.bam   1468223  0.535223
1-E06-F4_S60_L002_R1_001.fastqGene.bam_alignment.stats  3376137 1-E06-F4_S60_L002_R1_001.fastqGene_sorted.bam   1543492  0.457177
1-E07-F5_S61_L002_R1_001.fastqGene.bam_alignment.stats  518993  1-E07-F5_S61_L002_R1_001.fastqGene_sorted.bam   187987   0.362215
1-E08-F6_S62_L002_R1_001.fastqGene.bam_alignment.stats  5003471 1-E08-F6_S62_L002_R1_001.fastqGene_sorted.bam   2535015  0.506651
1-E09-F7_S63_L002_R1_001.fastqGene.bam_alignment.stats  4478648 1-E09-F7_S63_L002_R1_001.fastqGene_sorted.bam   2267976  0.506397
1-E10-F8_S64_L002_R1_001.fastqGene.bam_alignment.stats  4287261 1-E10-F8_S64_L002_R1_001.fastqGene_sorted.bam   2148095  0.501041
1-E11-F9_S65_L002_R1_001.fastqGene.bam_alignment.stats  4597089 1-E11-F9_S65_L002_R1_001.fastqGene_sorted.bam   2361679  0.513734
1-E12-F10_S66_L002_R1_001.fastqGene.bam_alignment.stats 5204689 1-E12-F10_S66_L002_R1_001.fastqGene_sorted.bam  2449930  0.470716
```
Hmm.  ~55% assignment is the highest, averaging around 50%.

#### How many of the total unique and multiple mapping reads were assigned to a gene feature?
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="UNPAIRED" {print $2}' )  <(cat *GeneMult_counts_genes.txt.summary|awk '$1=="Status" {print $2}' ) <(cat *GeneMult_counts_genes.txt.summary | awk '$1=="Assigned" {print $2}' ) <(cat *GeneMult_counts_genes.txt.summary | awk '$1=="Unassigned" {print $2}' ) |awk '{print $0,$4/$2}' |less

1-A01-A1_S7_L002_R1_001.fastqGene.bam_alignment.stats   6103950 1-A01-A1_S7_L002_R1_001.fastqGene_sorted.bam    3177782  0.520611
1-A02-A2_S8_L002_R1_001.fastqGene.bam_alignment.stats   3800050 1-A02-A2_S8_L002_R1_001.fastqGene_sorted.bam    2090815  0.550207
1-A03-A3_S9_L002_R1_001.fastqGene.bam_alignment.stats   3724472 1-A03-A3_S9_L002_R1_001.fastqGene_sorted.bam    2049499  0.550279
1-A04-A4_S10_L002_R1_001.fastqGene.bam_alignment.stats  4604921 1-A04-A4_S10_L002_R1_001.fastqGene_sorted.bam   2282808  0.495732
1-A05-A5_S11_L002_R1_001.fastqGene.bam_alignment.stats  3510078 1-A05-A5_S11_L002_R1_001.fastqGene_sorted.bam   1868040  0.532193
1-A06-A6_S12_L002_R1_001.fastqGene.bam_alignment.stats  4919113 1-A06-A6_S12_L002_R1_001.fastqGene_sorted.bam   2612014  0.530993
1-A07-A7_S13_L002_R1_001.fastqGene.bam_alignment.stats  5920091 1-A07-A7_S13_L002_R1_001.fastqGene_sorted.bam   3198690  0.540311
1-A08-A8_S14_L002_R1_001.fastqGene.bam_alignment.stats  6352631 1-A08-A8_S14_L002_R1_001.fastqGene_sorted.bam   3270028  0.514752
1-A09-A9_S15_L002_R1_001.fastqGene.bam_alignment.stats  6078054 1-A09-A9_S15_L002_R1_001.fastqGene_sorted.bam   3418230  0.562389
1-A10-A10_S16_L002_R1_001.fastqGene.bam_alignment.stats 5341438 1-A10-A10_S16_L002_R1_001.fastqGene_sorted.bam  2927251  0.548027
1-A11-A11_S17_L002_R1_001.fastqGene.bam_alignment.stats 5527935 1-A11-A11_S17_L002_R1_001.fastqGene_sorted.bam  2915608  0.527432
1-A12-A12_S18_L002_R1_001.fastqGene.bam_alignment.stats 5515285 1-A12-A12_S18_L002_R1_001.fastqGene_sorted.bam  2769828  0.502209
1-B01-A13_S19_L002_R1_001.fastqGene.bam_alignment.stats 6509294 1-B01-A13_S19_L002_R1_001.fastqGene_sorted.bam  3342529  0.513501
1-B02-A14_S20_L002_R1_001.fastqGene.bam_alignment.stats 3897336 1-B02-A14_S20_L002_R1_001.fastqGene_sorted.bam  1923324  0.493497
1-B03-A16_S21_L002_R1_001-024.fastqGene.bam_alignment.stats     23324186        1-B03-A16_S21_L002_R1_001-024.fastqGene_sorted.bam      11205040         0.480404
1-B04-A17_S22_L002_R1_001.fastqGene.bam_alignment.stats 3598206 1-B04-A17_S22_L002_R1_001.fastqGene_sorted.bam  1764183  0.490295
1-B05-A18_S23_L002_R1_001.fastqGene.bam_alignment.stats 2951360 1-B05-A18_S23_L002_R1_001.fastqGene_sorted.bam  1452629  0.49219
1-B06-A19_S24_L002_R1_001.fastqGene.bam_alignment.stats 5570408 1-B06-A19_S24_L002_R1_001.fastqGene_sorted.bam  2927654  0.525573
1-B07-B2_S25_L002_R1_001.fastqGene.bam_alignment.stats  6208488 1-B07-B2_S25_L002_R1_001.fastqGene_sorted.bam   3245555  0.522761
1-B08-B3_S26_L002_R1_001.fastqGene.bam_alignment.stats  5910328 1-B08-B3_S26_L002_R1_001.fastqGene_sorted.bam   2944766  0.498241
1-B09-B4_S27_L002_R1_001.fastqGene.bam_alignment.stats  4150147 1-B09-B4_S27_L002_R1_001.fastqGene_sorted.bam   2232897  0.538028
1-B10-B5_S28_L002_R1_001.fastqGene.bam_alignment.stats  5304293 1-B10-B5_S28_L002_R1_001.fastqGene_sorted.bam   2641090  0.497916
1-B11-B6_S29_L002_R1_001.fastqGene.bam_alignment.stats  2544152 1-B11-B6_S29_L002_R1_001.fastqGene_sorted.bam   930609   0.365784
1-B12-B7_S30_L002_R1_001.fastqGene.bam_alignment.stats  7076685 1-B12-B7_S30_L002_R1_001.fastqGene_sorted.bam   3738358  0.528264
1-C01-B8_S31_L002_R1_001.fastqGene.bam_alignment.stats  5274379 1-C01-B8_S31_L002_R1_001.fastqGene_sorted.bam   2957617  0.560752
1-C02-B10_S32_L002_R1_001.fastqGene.bam_alignment.stats 4254713 1-C02-B10_S32_L002_R1_001.fastqGene_sorted.bam  2164537  0.508739
1-C03-C1_S33_L002_R1_001.fastqGene.bam_alignment.stats  4210123 1-C03-C1_S33_L002_R1_001.fastqGene_sorted.bam   1946631  0.462369
1-C04-C2_S34_L002_R1_001.fastqGene.bam_alignment.stats  2958147 1-C04-C2_S34_L002_R1_001.fastqGene_sorted.bam   1503081  0.508116
1-C05-C3_S35_L002_R1_001.fastqGene.bam_alignment.stats  4199024 1-C05-C3_S35_L002_R1_001.fastqGene_sorted.bam   2325777  0.553885
1-C06-C4_S36_L002_R1_001.fastqGene.bam_alignment.stats  4783587 1-C06-C4_S36_L002_R1_001.fastqGene_sorted.bam   2677025  0.559627
1-C07-C6_S37_L002_R1_001.fastqGene.bam_alignment.stats  4191393 1-C07-C6_S37_L002_R1_001.fastqGene_sorted.bam   1981814  0.472829
1-C08-C7_S38_L002_R1_001.fastqGene.bam_alignment.stats  4940838 1-C08-C7_S38_L002_R1_001.fastqGene_sorted.bam   2449032  0.495671
1-C09-C8_S39_L002_R1_001.fastqGene.bam_alignment.stats  4717779 1-C09-C8_S39_L002_R1_001.fastqGene_sorted.bam   2442647  0.517754
1-C10-C9_S40_L002_R1_001.fastqGene.bam_alignment.stats  5864332 1-C10-C9_S40_L002_R1_001.fastqGene_sorted.bam   2922041  0.498273
1-C11-D1_S41_L002_R1_001.fastqGene.bam_alignment.stats  5446043 1-C11-D1_S41_L002_R1_001.fastqGene_sorted.bam   2929565  0.537925
1-C12-D2_S42_L002_R1_001.fastqGene.bam_alignment.stats  4443892 1-C12-D2_S42_L002_R1_001.fastqGene_sorted.bam   2329573  0.524219
1-D01-D5_S43_L002_R1_001.fastqGene.bam_alignment.stats  5688965 1-D01-D5_S43_L002_R1_001.fastqGene_sorted.bam   2685627  0.472077
1-D02-D6_S44_L002_R1_001.fastqGene.bam_alignment.stats  3583246 1-D02-D6_S44_L002_R1_001.fastqGene_sorted.bam   1853558  0.517285
1-D03-D7_S45_L002_R1_001.fastqGene.bam_alignment.stats  4803348 1-D03-D7_S45_L002_R1_001.fastqGene_sorted.bam   2027358  0.422072
1-D04-D8_S46_L002_R1_001.fastqGene.bam_alignment.stats  3265963 1-D04-D8_S46_L002_R1_001.fastqGene_sorted.bam   1368875  0.419134
1-D05-D9_S47_L002_R1_001.fastqGene.bam_alignment.stats  4167619 1-D05-D9_S47_L002_R1_001.fastqGene_sorted.bam   2025694  0.486055
1-D06-D10_S48_L002_R1_001.fastqGene.bam_alignment.stats 4350687 1-D06-D10_S48_L002_R1_001.fastqGene_sorted.bam  1944658  0.446977
1-D07-D11_S49_L002_R1_001.fastqGene.bam_alignment.stats 8144167 1-D07-D11_S49_L002_R1_001.fastqGene_sorted.bam  4219599  0.518113
1-D08-E1_S50_L002_R1_001.fastqGene.bam_alignment.stats  4785525 1-D08-E1_S50_L002_R1_001.fastqGene_sorted.bam   2110222  0.440959
1-D09-E2_S51_L002_R1_001.fastqGene.bam_alignment.stats  4628678 1-D09-E2_S51_L002_R1_001.fastqGene_sorted.bam   2342899  0.50617
1-D10-E5_S52_L002_R1_001.fastqGene.bam_alignment.stats  4738622 1-D10-E5_S52_L002_R1_001.fastqGene_sorted.bam   2104890  0.444199
1-D11-E7_S53_L002_R1_001.fastqGene.bam_alignment.stats  2707274 1-D11-E7_S53_L002_R1_001.fastqGene_sorted.bam   1454999  0.537441
1-D12-E8_S54_L002_R1_001.fastqGene.bam_alignment.stats  5761711 1-D12-E8_S54_L002_R1_001.fastqGene_sorted.bam   3203277  0.555959
1-E01-E9_S55_L002_R1_001.fastqGene.bam_alignment.stats  3718478 1-E01-E9_S55_L002_R1_001.fastqGene_sorted.bam   2027169  0.545161
1-E02-E10_S56_L002_R1_001.fastqGene.bam_alignment.stats 5871181 1-E02-E10_S56_L002_R1_001.fastqGene_sorted.bam  2807353  0.478158
1-E03-F1_S57_L002_R1_001.fastqGene.bam_alignment.stats  5837908 1-E03-F1_S57_L002_R1_001.fastqGene_sorted.bam   3197365  0.54769
1-E04-F2_S58_L002_R1_001.fastqGene.bam_alignment.stats  10028269        1-E04-F2_S58_L002_R1_001.fastqGene_sorted.bam   4841166  0.482752
1-E05-F3_S59_L002_R1_001.fastqGene.bam_alignment.stats  2743201 1-E05-F3_S59_L002_R1_001.fastqGene_sorted.bam   1494017  0.544625
1-E06-F4_S60_L002_R1_001.fastqGene.bam_alignment.stats  3376137 1-E06-F4_S60_L002_R1_001.fastqGene_sorted.bam   1572996  0.465916
1-E07-F5_S61_L002_R1_001.fastqGene.bam_alignment.stats  518993  1-E07-F5_S61_L002_R1_001.fastqGene_sorted.bam   191690   0.36935
1-E08-F6_S62_L002_R1_001.fastqGene.bam_alignment.stats  5003471 1-E08-F6_S62_L002_R1_001.fastqGene_sorted.bam   2580247  0.515691
1-E09-F7_S63_L002_R1_001.fastqGene.bam_alignment.stats  4478648 1-E09-F7_S63_L002_R1_001.fastqGene_sorted.bam   2309923  0.515763
1-E10-F8_S64_L002_R1_001.fastqGene.bam_alignment.stats  4287261 1-E10-F8_S64_L002_R1_001.fastqGene_sorted.bam   2186380  0.509971
1-E11-F9_S65_L002_R1_001.fastqGene.bam_alignment.stats  4597089 1-E11-F9_S65_L002_R1_001.fastqGene_sorted.bam   2403804  0.522897
1-E12-F10_S66_L002_R1_001.fastqGene.bam_alignment.stats 5204689 1-E12-F10_S66_L002_R1_001.fastqGene_sorted.bam  2496585  0.47968
```
Assigned ~37-56% of reads to a feature, including multi mapping reads.

#### How many of the total unique  mapping reads were assigned to a mRNA feature?
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="UNPAIRED" {print $2}' )  <(cat *mRNA_counts_genes.txt.summary|awk '$1=="Status" {print $2}' ) <(cat *mRNA_counts_genes.txt.summary | awk '$1=="Assigned" {print $2}' ) <(cat *mRNA_counts_genes.txt.summary | awk '$1=="Unassigned" {print $2}' ) |awk '{print $0,$4/$2}' |less
1-A01-A1_S7_L002_R1_001.fastqGene.bam_alignment.stats   6103950 1-A01-A1_S7_L002_R1_001.fastqGene_sorted.bam    1318529  0.216012
1-A02-A2_S8_L002_R1_001.fastqGene.bam_alignment.stats   3800050 1-A02-A2_S8_L002_R1_001.fastqGene_sorted.bam    916331   0.241137
1-A03-A3_S9_L002_R1_001.fastqGene.bam_alignment.stats   3724472 1-A03-A3_S9_L002_R1_001.fastqGene_sorted.bam    877356   0.235565
1-A04-A4_S10_L002_R1_001.fastqGene.bam_alignment.stats  4604921 1-A04-A4_S10_L002_R1_001.fastqGene_sorted.bam   999795   0.217114
1-A05-A5_S11_L002_R1_001.fastqGene.bam_alignment.stats  3510078 1-A05-A5_S11_L002_R1_001.fastqGene_sorted.bam   798419   0.227465
1-A06-A6_S12_L002_R1_001.fastqGene.bam_alignment.stats  4919113 1-A06-A6_S12_L002_R1_001.fastqGene_sorted.bam   1136302  0.230997
1-A07-A7_S13_L002_R1_001.fastqGene.bam_alignment.stats  5920091 1-A07-A7_S13_L002_R1_001.fastqGene_sorted.bam   1355258  0.228925
1-A08-A8_S14_L002_R1_001.fastqGene.bam_alignment.stats  6352631 1-A08-A8_S14_L002_R1_001.fastqGene_sorted.bam   1361569  0.214332
1-A09-A9_S15_L002_R1_001.fastqGene.bam_alignment.stats  6078054 1-A09-A9_S15_L002_R1_001.fastqGene_sorted.bam   1530919  0.251877
1-A10-A10_S16_L002_R1_001.fastqGene.bam_alignment.stats 5341438 1-A10-A10_S16_L002_R1_001.fastqGene_sorted.bam  1295229  0.242487
1-A11-A11_S17_L002_R1_001.fastqGene.bam_alignment.stats 5527935 1-A11-A11_S17_L002_R1_001.fastqGene_sorted.bam  1239846  0.224287
1-A12-A12_S18_L002_R1_001.fastqGene.bam_alignment.stats 5515285 1-A12-A12_S18_L002_R1_001.fastqGene_sorted.bam  1130230  0.204927
1-B01-A13_S19_L002_R1_001.fastqGene.bam_alignment.stats 6509294 1-B01-A13_S19_L002_R1_001.fastqGene_sorted.bam  1451153  0.222936
1-B02-A14_S20_L002_R1_001.fastqGene.bam_alignment.stats 3897336 1-B02-A14_S20_L002_R1_001.fastqGene_sorted.bam  853687   0.219044
1-B03-A16_S21_L002_R1_001-024.fastqGene.bam_alignment.stats     23324186        1-B03-A16_S21_L002_R1_001-024.fastqGene_sorted.bam      5329396  0.228492
1-B04-A17_S22_L002_R1_001.fastqGene.bam_alignment.stats 3598206 1-B04-A17_S22_L002_R1_001.fastqGene_sorted.bam  823575   0.228885
1-B05-A18_S23_L002_R1_001.fastqGene.bam_alignment.stats 2951360 1-B05-A18_S23_L002_R1_001.fastqGene_sorted.bam  677497   0.229554
1-B06-A19_S24_L002_R1_001.fastqGene.bam_alignment.stats 5570408 1-B06-A19_S24_L002_R1_001.fastqGene_sorted.bam  1330999  0.238941
1-B07-B2_S25_L002_R1_001.fastqGene.bam_alignment.stats  6208488 1-B07-B2_S25_L002_R1_001.fastqGene_sorted.bam   1468925  0.236599
1-B08-B3_S26_L002_R1_001.fastqGene.bam_alignment.stats  5910328 1-B08-B3_S26_L002_R1_001.fastqGene_sorted.bam   1254969  0.212335
1-B09-B4_S27_L002_R1_001.fastqGene.bam_alignment.stats  4150147 1-B09-B4_S27_L002_R1_001.fastqGene_sorted.bam   1004645  0.242075
1-B10-B5_S28_L002_R1_001.fastqGene.bam_alignment.stats  5304293 1-B10-B5_S28_L002_R1_001.fastqGene_sorted.bam   1121411  0.211416
1-B11-B6_S29_L002_R1_001.fastqGene.bam_alignment.stats  2544152 1-B11-B6_S29_L002_R1_001.fastqGene_sorted.bam   364141   0.143129
1-B12-B7_S30_L002_R1_001.fastqGene.bam_alignment.stats  7076685 1-B12-B7_S30_L002_R1_001.fastqGene_sorted.bam   1450912  0.205027
1-C01-B8_S31_L002_R1_001.fastqGene.bam_alignment.stats  5274379 1-C01-B8_S31_L002_R1_001.fastqGene_sorted.bam   1162905  0.220482
1-C02-B10_S32_L002_R1_001.fastqGene.bam_alignment.stats 4254713 1-C02-B10_S32_L002_R1_001.fastqGene_sorted.bam  846536   0.198964
1-C03-C1_S33_L002_R1_001.fastqGene.bam_alignment.stats  4210123 1-C03-C1_S33_L002_R1_001.fastqGene_sorted.bam   801180   0.190298
1-C04-C2_S34_L002_R1_001.fastqGene.bam_alignment.stats  2958147 1-C04-C2_S34_L002_R1_001.fastqGene_sorted.bam   640381   0.21648
1-C05-C3_S35_L002_R1_001.fastqGene.bam_alignment.stats  4199024 1-C05-C3_S35_L002_R1_001.fastqGene_sorted.bam   955275   0.227499
1-C06-C4_S36_L002_R1_001.fastqGene.bam_alignment.stats  4783587 1-C06-C4_S36_L002_R1_001.fastqGene_sorted.bam   1102907  0.230561
1-C07-C6_S37_L002_R1_001.fastqGene.bam_alignment.stats  4191393 1-C07-C6_S37_L002_R1_001.fastqGene_sorted.bam   892108   0.212843
1-C08-C7_S38_L002_R1_001.fastqGene.bam_alignment.stats  4940838 1-C08-C7_S38_L002_R1_001.fastqGene_sorted.bam   994842   0.201351
1-C09-C8_S39_L002_R1_001.fastqGene.bam_alignment.stats  4717779 1-C09-C8_S39_L002_R1_001.fastqGene_sorted.bam   1086023  0.230198
1-C10-C9_S40_L002_R1_001.fastqGene.bam_alignment.stats  5864332 1-C10-C9_S40_L002_R1_001.fastqGene_sorted.bam   1219582  0.207966
1-C11-D1_S41_L002_R1_001.fastqGene.bam_alignment.stats  5446043 1-C11-D1_S41_L002_R1_001.fastqGene_sorted.bam   1266411  0.232538
1-C12-D2_S42_L002_R1_001.fastqGene.bam_alignment.stats  4443892 1-C12-D2_S42_L002_R1_001.fastqGene_sorted.bam   960213   0.216075
1-D01-D5_S43_L002_R1_001.fastqGene.bam_alignment.stats  5688965 1-D01-D5_S43_L002_R1_001.fastqGene_sorted.bam   1126529  0.19802
1-D02-D6_S44_L002_R1_001.fastqGene.bam_alignment.stats  3583246 1-D02-D6_S44_L002_R1_001.fastqGene_sorted.bam   784411   0.218911
1-D03-D7_S45_L002_R1_001.fastqGene.bam_alignment.stats  4803348 1-D03-D7_S45_L002_R1_001.fastqGene_sorted.bam   874329   0.182025
1-D04-D8_S46_L002_R1_001.fastqGene.bam_alignment.stats  3265963 1-D04-D8_S46_L002_R1_001.fastqGene_sorted.bam   557427   0.170678
1-D05-D9_S47_L002_R1_001.fastqGene.bam_alignment.stats  4167619 1-D05-D9_S47_L002_R1_001.fastqGene_sorted.bam   892148   0.214067
1-D06-D10_S48_L002_R1_001.fastqGene.bam_alignment.stats 4350687 1-D06-D10_S48_L002_R1_001.fastqGene_sorted.bam  823976   0.18939
1-D07-D11_S49_L002_R1_001.fastqGene.bam_alignment.stats 8144167 1-D07-D11_S49_L002_R1_001.fastqGene_sorted.bam  1845484  0.226602
1-D08-E1_S50_L002_R1_001.fastqGene.bam_alignment.stats  4785525 1-D08-E1_S50_L002_R1_001.fastqGene_sorted.bam   878859   0.183649
1-D09-E2_S51_L002_R1_001.fastqGene.bam_alignment.stats  4628678 1-D09-E2_S51_L002_R1_001.fastqGene_sorted.bam   1033903  0.223369
1-D10-E5_S52_L002_R1_001.fastqGene.bam_alignment.stats  4738622 1-D10-E5_S52_L002_R1_001.fastqGene_sorted.bam   917049   0.193527
1-D11-E7_S53_L002_R1_001.fastqGene.bam_alignment.stats  2707274 1-D11-E7_S53_L002_R1_001.fastqGene_sorted.bam   609019   0.224957
1-D12-E8_S54_L002_R1_001.fastqGene.bam_alignment.stats  5761711 1-D12-E8_S54_L002_R1_001.fastqGene_sorted.bam   1301351  0.225862
1-E01-E9_S55_L002_R1_001.fastqGene.bam_alignment.stats  3718478 1-E01-E9_S55_L002_R1_001.fastqGene_sorted.bam   859492   0.231141
1-E02-E10_S56_L002_R1_001.fastqGene.bam_alignment.stats 5871181 1-E02-E10_S56_L002_R1_001.fastqGene_sorted.bam  1137887  0.193809
1-E03-F1_S57_L002_R1_001.fastqGene.bam_alignment.stats  5837908 1-E03-F1_S57_L002_R1_001.fastqGene_sorted.bam   1353179  0.231792
1-E04-F2_S58_L002_R1_001.fastqGene.bam_alignment.stats  10028269        1-E04-F2_S58_L002_R1_001.fastqGene_sorted.bam   2041581  0.203583
1-E05-F3_S59_L002_R1_001.fastqGene.bam_alignment.stats  2743201 1-E05-F3_S59_L002_R1_001.fastqGene_sorted.bam   670285   0.244344
1-E06-F4_S60_L002_R1_001.fastqGene.bam_alignment.stats  3376137 1-E06-F4_S60_L002_R1_001.fastqGene_sorted.bam   656625   0.19449
1-E07-F5_S61_L002_R1_001.fastqGene.bam_alignment.stats  518993  1-E07-F5_S61_L002_R1_001.fastqGene_sorted.bam   80541    0.155187
1-E08-F6_S62_L002_R1_001.fastqGene.bam_alignment.stats  5003471 1-E08-F6_S62_L002_R1_001.fastqGene_sorted.bam   1110103  0.221867
1-E09-F7_S63_L002_R1_001.fastqGene.bam_alignment.stats  4478648 1-E09-F7_S63_L002_R1_001.fastqGene_sorted.bam   979317   0.218664
1-E10-F8_S64_L002_R1_001.fastqGene.bam_alignment.stats  4287261 1-E10-F8_S64_L002_R1_001.fastqGene_sorted.bam   941552   0.219616
1-E11-F9_S65_L002_R1_001.fastqGene.bam_alignment.stats  4597089 1-E11-F9_S65_L002_R1_001.fastqGene_sorted.bam   1076598  0.234191
1-E12-F10_S66_L002_R1_001.fastqGene.bam_alignment.stats 5204689 1-E12-F10_S66_L002_R1_001.fastqGene_sorted.bam  1056320  0.202955
```
Only unique reads assigned at about 19-25% of mRNA features.
