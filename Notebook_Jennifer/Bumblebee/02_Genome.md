# 02 Genome

Fetch bumblebee `Bombus impatiens` reference files from NCBI

```
mkdir 02_Genome; cd 02_Genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_rna.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_genomic.gff.gz
```

Notes: Should be using reference files from the hymenoptera site

* https://hymenoptera.elsiklab.missouri.edu/genome_fasta

but where are the gff and transcriptome files?

## Pull out the Uniprot IDS

```
cat GCF_000188095.3_BIMP_2.2_genomic.gff |awk -F'\t' '$3=="gene" {print $9}' |awk -F';' '{print $1,$2}'  > loc_to_uniprot.txt
```

```
less loc_to_uniprot.txt

ID=gene-LOC100740276 Dbxref=GeneID:100740276
ID=gene-LOC100740157 Dbxref=GeneID:100740157
ID=gene-LOC100742884 Dbxref=GeneID:100742884
ID=gene-LOC100740399 Dbxref=GeneID:100740399
ID=gene-LOC100740519 Dbxref=GeneID:100740519
ID=gene-LOC100743001 Dbxref=GeneID:100743001
ID=gene-LOC100740639 Dbxref=GeneID:100740639
....
```

Might need to switch to the hymenoptera site, found the gff

```
wget https://hymenoptera.elsiklab.missouri.edu/sites/hymenoptera.org/files/data/gff3/ogs/Bombus_impatiens_bimp_OGSv1.0_liftover_BIMP_2.2.gff3

less Bombus_impatiens_bimp_OGSv1.0_liftover_BIMP_2.2.gff3
```

Which list `BIMP` IDs

```
NT_179809.1	bimp_OGSv1.0	gene	94	591	1	-	.	ID=BIMP10001;Name=BIMP10001;score=1;source_info=Official Gene Set;
NT_179809.1	bimp_OGSv1.0	mRNA	94	591	1	-	.	ID=BIMP10001-RA;Name=BIMP10001-RA;Parent=BIMP10001;score=1;source_info=Official Gene Set;
NT_179809.1	bimp_OGSv1.0	exon	94	591	.	-	0	ID=exon1;Parent=BIMP10001-RA;source_info=Official Gene Set;
NT_179809.1	bimp_OGSv1.0	CDS	94	591	.	-	0	ID=cds1;Parent=BIMP10001-RA;source_info=Official Gene Set;
NT_178602.1	bimp_OGSv1.0	gene	651	764	0.999936	+	.	ID=BIMP10002;Name=BIMP10002;score=0.999936;source_info=Official Gene Set;
NT_178602.1	bimp_OGSv1.0	mRNA	651	764	0.999936	+	.	ID=BIMP10002-RA;Name=BIMP10002-RA;Parent=BIMP10002;score=0.999936;source_info=Official Gene Set;
...
```
