# Install and run rseqc to determine stranding information for maize paired end.

```
#/work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align

ml miniconda3
conda create -n rseqc

source activate rseqc

#this did not finish, using pip
conda install -c bioconda rseq

pip install --user RSeQC




```


### Convert maize gff to bed12 format for rseqc
```
ml bedops
gff2bed <GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff >maizeGFF.bed12
```

### Identify stranding information
```
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.5082
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2406
Fraction of reads explained by "1+-,1-+,2++,2--": 0.2512
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573505_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3724
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3053
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3222
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573506_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3975
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2955
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3070
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573507_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3976
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2977
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3047
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573508_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3600
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3130
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3269
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573509_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.5146
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2338
Fraction of reads explained by "1+-,1-+,2++,2--": 0.2516
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573510_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled

This is PairEnd Data
Fraction of reads failed to determine: 0.3387
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3272
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3341
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573511_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3559
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3160
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3281
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573512_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.4640
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2660
Fraction of reads explained by "1+-,1-+,2++,2--": 0.2700
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573513_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3429
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3106
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3464
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573514_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.4062
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2726
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3211
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573515_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3973
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2846
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3180
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573516_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3291
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3121
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3588
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573517_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3857
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2775
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3368
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573518_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.4687
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2493
Fraction of reads explained by "1+-,1-+,2++,2--": 0.2819
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573519_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled

This is PairEnd Data
Fraction of reads failed to determine: 0.2921
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3463
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3617
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573520_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3442
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3242
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3315
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573521_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3660
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3126
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3214
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573522_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3827
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3056
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3117
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573523_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.4076
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2925
Fraction of reads explained by "1+-,1-+,2++,2--": 0.2999
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573524_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3597
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3141
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3261
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573525_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3001
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3460
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3540
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573526_1Gene_sorted.bam'

Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3692
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3142
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3166
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573527_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3451
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3267
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3282

```
