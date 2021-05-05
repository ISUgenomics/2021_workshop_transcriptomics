# 11 Map to Gene Names

Pull out `rna_list.txt` from the rRNA file

```
zcat GCF_000188095.3_BIMP_2.2_rna.fna.gz | grep "^>" > rna_list.txt

less rna_list.txt

>NM_001280122.1 Bombus impatiens cryptochrome 2 (CRY2), mRNA
>NM_001280123.1 Bombus impatiens ultraviolet-sensitive opsin (SWRh), mRNA
>NM_001365858.1 Bombus impatiens cytochrome P450 6a14-like (LOC100745439), mRNA
>NM_001365859.1 Bombus impatiens cytochrome P450 6a14 (LOC113218476), mRNA
>XM_003484340.4 PREDICTED: Bombus impatiens PHD finger-like domain-containing protein 5A (LOC100740157), mRNA
```

<!-- NCBI option, commented out because I didn't use this section
Download gene_result.txt from NCBI

* https://www.ncbi.nlm.nih.gov/gene/?term=txid132113[Organism:noexp]


```
less gene_result.txt

tax_id Org_name GeneID CurrentID Status Symbol Aliases description other_designations map_location chromosome genomic_nucleotide_accession.version start_position_on_the_genomic_accession end_position_on_the_genomic_accession orientation exon_count OMIM 
132113 Bombus impatiens 100607971 0 live SWRh  ultraviolet-sensitive opsin ultraviolet-sensitive opsin  Un       
132113 Bombus impatiens 100607970 0 live CRY2  cryptochrome 2 cryptochrome 2  Un       
132113 Bombus impatiens 100743205 0 live Osi6  Osiris 6 uncharacterized protein LOC100743205  Un       
132113 Bombus impatiens 100743092 0 live LOC100743092 BimIRP30, IRP30 slit homolog 3 protein slit homolog 3 protein  Un     
...
```

Hmm... `gene_result.txt` only has 13171 rows, while `rna_list.txt` has 27810 rows. Maybe continue with `rna_list.txt`.

-->

<details><summary>reformat.pl</summary>

```
#! /usr/bin/env perl

use strict;
use warnings;

while(<>){
    chomp;
    if(/>(\S+)\s.+impatiens (\S.+)\((\S+)\)(.+)/){
        print("$1\t$3\t$2($3)$4\n")
    }else{
        print("$_\n")
    }
}
```

</details>

```
./reformat.pl rna_ids.txt > id_gene.txt

less id_gene.txt

NM_001280122.1	CRY2	cryptochrome 2 (CRY2), mRNA
NM_001280123.1	SWRh	ultraviolet-sensitive opsin (SWRh), mRNA
NM_001365858.1	LOC100745439	cytochrome P450 6a14-like (LOC100745439), mRNA
NM_001365859.1	LOC113218476	cytochrome P450 6a14 (LOC113218476), mRNA
XM_003484340.4	LOC100740157	PHD finger-like domain-containing protein 5A (LOC100740157), mRNA
...
```

Merge with Kallisto and Salmon results

```
library(tidyverse)
library(magrittr)

data <- readxl::read_excel("salmon_counts.xlsx", sheet = "NumReads")
head(data)

geneids <- read_delim("../id_gene.txt", delim="\t", col_names = c("id", "gene", "description"))
row.names(geneids) = geneids$id

sample_names <- names(data)[-1]

data$gene = geneids[data$Name,]$gene
data$description = geneids[data$Name,]$description

cdata<-data %>%
  select(Name, gene, description, all_of(sample_names))

write_delim(cdata, "salmon_NumReads.txt", delim="\t")
```

Merge with GSNAP DESeq2 results

```
library(tidyverse)
library(magrittr)

data <- readxl::read_excel("gsnap_DESeq2_results.xlsx", sheet = "mRNA_DESeq2")

geneids <- read_delim("../id_gene.txt", delim="\t", col_names = c("mrna","gene", "description"))

row.names(geneids) = paste("rna-",geneids$mrna,sep="")

sample_names <- names(data)[-1]

#data$gene = geneids[data$Gene,]$gene
data$description = geneids[data$Gene,]$description

cdata<-data %>%
  select(Gene, description, all_of(sample_names))

write_delim(cdata, "gsnap_mRNA_DESeq2.txt", delim="\t")
```


