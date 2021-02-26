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

Although

```


#### Genome and Annotation
`/work/gif/Siva/Transcriptomic_Worshop/Genome`

```
cd /work/gif/Siva/Transcriptomic_Worshop/Genome/
mkdir Ver4_release50
cd Ver4_release50/

#Genome:
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
# annotation
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-50/gtf/zea_mays/Zea_mays.B73_RefGen_v4.50.gtf.gz

mkdir Ver5
cd Ver5
#Genome:
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
#Annotation:
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gene.fa.gz

 ```

 #### Hisat2 index
* Built genome index for both ver 4 and ver5(NAM)
 ```
cd Ver4_release50/
gunzip *
module load hisat2
hisat2-build -p18 Zea_mays.B73_RefGen_v4.dna.toplevel.fa Zea_mays_Ver4_50 >& ver4_build.log&

 cd ../Ver5/
gunzip *
hisat2-build -p 18 Zm-B73-REFERENCE-NAM-5.0.fa Zm-B73-REFERENCE-NAM-ver5 >& ver5_build.log&
 ```


* Version 4 Index details:
 ```
 Returning from HGFM constructor
Headers:
    len: 2104350182
    gbwtLen: 2104350183
    nodes: 2104350183
    sz: 526087546
    gbwtSz: 526087546
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 0
    eftabSz: 0
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 131521887
    offsSz: 526087548
    lineSz: 64
    sideSz: 64
    sideGbwtSz: 48
    sideGbwtLen: 192
    numSides: 10960158
    numLines: 10960158
    gbwtTotLen: 701450112
    gbwtTotSz: 701450112
    reverse: 0
    linearFM: Yes
Total time for call to driver() for forward index: 00:08:03
 ```
* Version 5 Index Details:
 ```
 Headers:
    len: 2178268108
    gbwtLen: 2178268109
    nodes: 2178268109
    sz: 544567027
    gbwtSz: 544567028
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 0
    eftabSz: 0
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 136141757
    offsSz: 544567028
    lineSz: 64
    sideSz: 64
    sideGbwtSz: 48
    sideGbwtLen: 192
    numSides: 11345147
    numLines: 11345147
    gbwtTotLen: 726089408
    gbwtTotSz: 726089408
    reverse: 0
    linearFM: Yes
Total time for call to driver() for forward index: 00:17:45
 ```
#### Hisat2 index
* With Ver 4 genome:

 ```
[csiva@nova 02_align]$ more hisat2.sh
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=36
#SBATCH --time=24:00:00
#SBATCH --job-name=Hisat2
#SBATCH --output=Hisat2.%j.out
#SBATCH --error=Hisat2.%j.err
#SBATCH --mail-user=csiva@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
set -o xtrace
cd ${SLURM_SUBMIT_DIR}
scontrol show job ${SLURM_JOB_ID}
# purge and load relevant modules.
module purge
module load parallel
module load hisat2
module load samtools
parallel -j4 "hisat2 -p 8 -x /work/gif/Siva/Transcriptomic_Worshop/Genome/Ver4_release50/Zea_mays_Ver4_50
 -1 {1} -2 {2} -S {1/.}.sam 2> {1/.}.log" ::: ../maize_data/*1.fastq :::+ ../maize_data/*2.fastq
parallel -j4 "samtools view --threads 8 -bS -o {/.}.bam {}" ::: *.sam
rm *.sam
```

#### Stringtie2 Assembly:
"StringTie is a fast and highly efficient assembler of RNAseq alignments into potential transcripts" ([StringTie Manual](https://ccb.jhu.edu/software/stringtie/index.shtml)). The version on Nova is old `StringTie v1.3.4a` so decided to download precompiled [binary](http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.4.Linux_x86_64.tar.gz) in `/work/gif/Siva/software` and `mv <binary> ~/.local/bin`. The version I now have is `StringTie v2.1.4`.

* Assembling the mapped reads to transcripts:

```
[csiva@nova008 03_Stringtie_Assembly]$ pwd
/work/gif/Siva/Transcriptomic_Worshop/03_Stringtie_Assembly
[csiva@nova008 03_Stringtie_Assembly]$ mkdir Assembly
```
* Sort the bam files:

```
/work/gif/Siva/Transcriptomic_Worshop/02_align/ver4
[csiva@nova008 ver4]$ ls

LOGS    SRR1573504_1.bam  SRR1573506_1.bam  SRR1573508_1.bam  SRR1573510_1.bam  SRR1573512_1.bam  SRR1573514_1.bam  SRR1573516_1.bam  SRR1573518_1.bam  SRR1573520_1.bam  SRR1573522_1.bam  SRR1573524_1.bam  SRR1573526_1.bam
sorted  SRR1573505_1.bam  SRR1573507_1.bam  SRR1573509_1.bam  SRR1573511_1.bam  SRR1573513_1.bam  SRR1573515_1.bam  SRR1573517_1.bam  SRR1573519_1.bam  SRR1573521_1.bam  SRR1573523_1.bam  SRR1573525_1.bam  SRR1573527_1.bam

parallel -j4 "samtools sort -@7 -o sorted/{.}.sorted.bam {}" ::: *bam
```
After sorting the bam files, the unsorted BAM files can be removed.

* Assembling known transcripts using known annotation file as a guide. `-G <GTF_file>`

```
parallel -j4 "stringtie {} -p8 -G ../Genome/Ver4_release50/Zea_mays.B73_RefGen_v4.50.gtf -o Assembly/{/.}.gtf 2> Assembly/{/.}.log" ::: ../02_align/ver4/sorted/*bam >&stringtie-run.log&
```
Files Location `/work/gif/Siva/Transcriptomic_Worshop/03_Stringtie_Assembly/Assembly`


* Install gffcompare

```
cd /work/gif/Siva/software/
git clone https://github.com/gpertea/gffcompare
cd gffcompare/
make
cp gffcompare ~/.local/bin/

```

* Merge all the stringtie transcripts:

```
ls *gtf | more > mergelist.txt

stringtie --merge -p 32 -G /work/gif/Siva/Transcriptomic_Worshop/Genome/Ver4_release50/Zea_mays.B73_RefGen_v4.50.gtf -o stringtie_merged_transcripts.gtf mergelist.txt >& stringtie-merge.log&

cat stringtie_merged_transcripts.gtf | head -4

# stringtie --merge -p 32 -G /work/gif/Siva/Transcriptomic_Worshop/Genome/Ver4_release50/Zea_mays.B73_RefGen_v4.50.gtf -o stringtie_merged_transcripts.gtf mergelist.txt
# StringTie version 2.1.4
1       StringTie       transcript      44289   49868   1000    +       .       gene_id "MSTRG.1"; transcript_id "Zm00001d027230_T001"; ref_gene_id "Zm00001d027230";
1       StringTie       exon    44289   44947   1000    +       .       gene_id "MSTRG.1"; transcript_id "Zm00001d027230_T001"; exon_number "1"; ref_gene_id "Zm00001d027230";
```
* Check the total number of transcripts assembled:
```
cat stringtie_merged_transcripts.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l
182912 ## Number of transcripts:
```

* Compare the stringtie assembly to known transcripts:
```
gffcompare -r /work/gif/Siva/Transcriptomic_Worshop/Genome/Ver4_release50/Zea_mays.B73_RefGen_v4.50.gtf -G -o merged stringtie_merged_transcripts.gtf
#gffcompare -r /work/gif/Siva/Transcriptomic_Worshop/Genome/Ver4_release50/Zea_mays.B73_RefGen_v4.50.gtf -G -o merged stringtie_merged_transcripts.gtf
#
```

* merged.stats
```
#= Summary for dataset: stringtie_merged_transcripts.gtf
#     Query mRNAs :  182912 in   53607 loci  (162112 multi-exon transcripts)
#            (20800 multi-transcript loci, ~3.4 transcripts per locus)
# Reference mRNAs :  134190 in   45116 loci  (120908 multi-exon)
# Super-loci w/ reference transcripts:    44581
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    89.3    |
        Exon level:    92.4     |    81.2    |
      Intron level:   100.0     |    87.3    |
Intron chain level:    99.9     |    74.5    |
  Transcript level:    99.5     |    73.0    |
       Locus level:    99.4     |    82.7    |

     Matching intron chains:  120730
       Matching transcripts:  133580
              Matching loci:   44835

          Missed exons:       0/403898  (  0.0%)
           Novel exons:   24796/449809  (  5.5%)
        Missed introns:       3/244908  (  0.0%)
         Novel introns:   15996/280444  (  5.7%)
           Missed loci:       0/45116   (  0.0%)
            Novel loci:    9026/53607   ( 16.8%)

 Total union super-loci across all input datasets: 53607
182912 out of 182912 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)
```

* Use sorted BAM files as input to perform transcript abundance estimation using the `-e` `-B` and `-G` options:

```
[csiva@nova008 04_Ballgown]$ pwd
/work/gif/Siva/Transcriptomic_Worshop/04_Ballgown
[csiva@nova008 04_Ballgown]$ parallel -j4 "stringtie -e -B -p 8 -o {/.}/{/.}.ball.gtf -G ../03_Stringtie_Assembly/Assembly/GTF/Merged_files/stringtie_merged_transcripts.gtf {} 2> {/.}.ball.log" ::: ../02_align/ver4/sorted/SRR15735*bam >& std.out&
```
Files snippet within each folder made for Ballgown`[csiva@nova001 04_Ballgown]$ pwd
/work/gif/Siva/Transcriptomic_Worshop/04_Ballgown`
```
[csiva@nova001 04_Ballgown]$ ls -1 SRR1573504_1.sorted/
e2t.ctab
e_data.ctab
i2t.ctab
i_data.ctab
SRR1573504_1.sorted.ball.gtf
t_data.ctab
```
There are a total of 24 folders within the `04_Ballgown` directory; total of 2G. I copied the entire directory to my local machine for working on differential expression using Ballgown R package.

See the [R script](/Notebook_Siva/Maize_ligule.R) for initial steps in Differential Expression analysis with this dataset.
