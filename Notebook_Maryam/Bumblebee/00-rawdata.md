# Raw data for BumbleBee

* Feb 2021
* Nova : /work/gif/Maryam/projects/Transcriptomics/00-rawdata/BumbleBee

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
* Note: I figured that the gff file I downloaded have the RefSeq annotations. So I have to use another gff file :

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/Gnomon_models/GCF_000188095.3_BIMP_2.2_gnomon_model.gff.gz
```
