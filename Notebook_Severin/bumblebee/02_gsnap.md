# run gsnap on Bumble Bee

* Feb 18, 2021
* /work/gif/TranscriptomicsWorkshop/severin/Bee/02_gsnap

## Softlink raw data and genome and gff

```
ls /work/gif/TranscriptomicsWorkshop/bumbleBee/*gz  | xargs -I xx ln -s xx
ln -s /work/gif/TranscriptomicsWorkshop/bumbleBee/01_Genome/GCF_000188095.3_BIMP_2.2_genomic.fna
cp /work/gif/TranscriptomicsWorkshop/bumbleBee/01_Genome/GCF_000188095.3_BIMP_2.2_genomic.gff.gz .
gunzip GCF_000188095.3_BIMP_2.2_genomic.gff.gz

```


### Build the Maize database for gsnap/gmap

* /work/gif/TranscriptomicsWorkshop/severin/Bee/02_gsnap

```
module load gmap-gsnap  # version 2019-05-12-zjqshxf
gmap_build -D BB -d BB GCF_000188095.3_BIMP_2.2_genomic.gff
```
