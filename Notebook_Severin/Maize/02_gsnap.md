# Align reads

* Feb 16, 2021
* /work/gif/TranscriptomicsWorkshop/severin/Maize/02_gsnap/

# Maize

### Softlink the files

```
ls /work/gif/TranscriptomicsWorkshop/Maize/*fastq | more | xargs -I xx ln -s xx
```

### Download the Maize genome

* /work/gif/TranscriptomicsWorkshop/severin/Maize/02_gsnap/assembly
Download and combine chromosomes

```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.*gz
zcat *.gz > Zea_All.fasta
```

### Build the Maize database for gsnap/gmap

*  /work/gif/TranscriptomicsWorkshop/severin/Maize/02_gsnap

```
gmap_build -D Zea -d B73 Zea_All.fasta
```
