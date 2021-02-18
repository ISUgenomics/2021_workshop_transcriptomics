# Alignment with Star Bee dataset
* Feb 2021
* Nova : /work/gif/Maryam/projects/Transcriptomics/02Bee-STAR

#### STAR aligner

[STAR](https://github.com/alexdobin/STAR) is an alignment tool designed to specifically address challenges of RNA-seq mapping using a strategy to account for spliced alignment. It has shown high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed.

* Cons:
Memory intensive
* Pros:
Fast mapping speed
* Strategy:
Tow-step mapping process:

1. Seed searching
2. Clustering, stitching, and scoring

### Ref data set:
* **Note**: Use the link bellow. [BeeBase](https://hymenoptera.elsiklab.missouri.edu/beebase) is been retired!

[Hymenoptera Genome Database](https://hymenoptera.elsiklab.missouri.edu/genome_fasta)
* Genome :
```
wget https://hymenoptera.elsiklab.missouri.edu/sites/hymenoptera.org/files/data/genomes/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa.gz
```
* gff file (we start with the protein coding gff3 file)
```
wget https://hymenoptera.elsiklab.missouri.edu/sites/hymenoptera.org/files/data/gff3/refseq/Bombus_impatiens_BIMP_2.2_RefSeq_proteincoding.gff3
```


#### First step: Generating genome index
Basic options:
```
--runThreadN  Number of threads
--runMode     genome Generate
--genomeDir   Path to the output genome indexes
--genomeFastaFiles  path to the genome fasta files
--sjdbGTFfile        path to annotations.gtf (Optional but highly recommended)
--sjdbOverhang      ReadLength-1
```

Building the index:
```bash
salloc -N1 -n8 -p DEBUG -t 1:00:00
cd /work/gif/Maryam/projects/Transcriptomics/02Bee-STAR
module load star

mv /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fasta

mkdir BeeDB

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./BeeDB  --sjdbGTFfile /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/Bombus_impatiens_BIMP_2.2_RefSeq_proteincoding.gff3 --genomeFastaFiles /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fasta --sjdbOverhang  49


```
Notes :
* I checked one read file and its length was 50 bp. So I choose 49 for the overhanging value int he options above!
* STAR does not build the output directory. it should be built before running STAR!
* I renamed the fasta file (changed `.fa` to `.fast`. I don't think it mattered.)

Error message:
`terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check
`

I found [this](https://groups.google.com/g/rna-star/c/ONinttkq1q0?pli=1):

```
it is possible to use the gff file, however you would need to specify --sjdbGTFtagExonParentTranscript parameter, which is typically 'Parent' for gff files -
that assigns exon to a transcript.
If you want to use --sjdbGTFtagExonParentTranscript GeneCounts option, the gff file would have to contain the gene attrbiute for each exon, specified in --sjdbGTFtagExonParentGene. Since typically this attribute is not present, it's best to convert the file to GTF with gene_id and transcript_id for each exons. You can use your own script, or gffread from Cufflinks package:
gffread -T small.gff3 -o small.gtf

```

Lets try again:
```bash
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./BeeDB  --sjdbGTFfile /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/Bombus_impatiens_BIMP_2.2_RefSeq_proteincoding.gff3 --genomeFastaFiles /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fasta --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang  49
```

It worked!

I have the genome indexed now!

### Second Step: Alignment

Script to run the STAR alignments :
* run_star.sub

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


echo >> $filename
echo "module load star"  >> $filename
echo >> $filename

echo "star --runThreadN 12 --readFilesCommand zcat --readFilesIn /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee/$read1  --genomeDir /work/gif/Maryam/projects/Transcriptomics/02Bee-STAR/BeeDB  --outFileNamePrefix StarOut" >> $filename
echo "scontrol show job \$SLURM_JOB_ID" >> $filename
```

running it within a for loop
```bash
for n in ../00-rawdata/BumbleBee/*.fastq.gz; do  ./run_star.sub $n ; done
```

* Notes:
1. Make sure you have the `log` directory when you run the script!
2. Alignments run so fast that I decided to run them all in an interactive node instead of submitting the jobs through slurm!

```bash
for n in ../00-rawdata/BumbleBee/*.fastq.gz; do   read1=$n; filename=`echo $read1 | cut -d'/' -f 4|  cut -d'_' -f 1` ; STAR --runThreadN 12  --readFilesIn $read1   --genomeDir /work/gif/Maryam/projects/Transcriptomics/02Bee-STAR/BeeDB  --readFilesCommand zcat  --outFileNamePrefix StarOut-$filename; done

```
