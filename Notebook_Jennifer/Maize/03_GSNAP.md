# 03 GSNAP

(1) Index reference

```
module load gmap-gsnap samtools subread

GENOME=02_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
G_NAME=$(basename $GENOME) |sed 's/.fna.gz//g'

gmap_build --gunzip -d ${G_NAME} -D gmapdb ${GENOME}
```

(2) Map reads

```
PROC=16

for READ in 00_Raw-Data/*1.fastq.gz
do
  READNAME=$(basename ${FILE}) |sed 's/.fastq.gz//g'
  READ2=$(echo $READ | sed 's/1.fastq.gz/2.fastq.gz/g')
  
  gsnap --gunzip -d ${G_NAME} -D gmapdb -N 1 -t 16 -B 4 -m 5 \
    --input-buffer-size=1000000 --output-buffer-size=1000000 \
	-A sam ${READ} ${READ2} | \
	samtools view --threads 16 -bS - > ${READNAME}.bam
  
done
```

(3) Quantify gene counts

```
GFF=02_Genome/GCF_000188095.3_BIMP_2.2_genomic.gff.gz

for BAM in 00_Raw-Data/*.bam
do
  READNAME=$(basename ${FILE}) |sed 's/.fastq.gz//g'
  
  featureCounts -T 16 -t gene -g ID -a ${GFF} -o ${READNAME}_genecounts.txt ${BAM}
  featureCounts -T 16 -t mRNA -g ID -a ${GFF} -o ${READNAME}_mRNAcounts.txt ${BAM}
  featureCounts -T 16 -M -t gene -g ID -a ${GFF} -o ${READNAME}_geneMultcounts.txt ${BAM}
done
```

(4) Combine gene counts into one file

```
#! /usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(readxl)

# === Get list of featureCount output files
featureCount_files <- list.files(path = "03_GSNAP", pattern = "*genecounts.txt$", full.names = TRUE)

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" ) %>%
  select(Geneid, ends_with(".bam")) %>%             
  pivot_longer(cols=ends_with(".bam")) %>% 
  mutate(
    name = gsub(".bam", "", name)
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(Geneid, ends_with(".bam")) %>%
    pivot_longer(cols=ends_with(".bam")) %>%
    mutate(
      name = gsub(".bam", "", name)
    )
  data = rbind(data, temp)
}

# === Convert to excel like data (wider)
wide_data <- data %>%
  pivot_wider(id_cols=Geneid)

# === Save tab delimited file (smaller file size)
write_delim(wide_data, "gsnap_genecounts.txt", delim="\t")

# === Save Excel file (can be easier to work with)
# writexl::write_xlsx(wide_data, "gsnap_genecounts.xlsx")
```

(5) Repeat above for `*mRNAcounts.txt` and `*geneMultcounts.txt` and combine all counts into one Excel file [`gsnap_counts.xlsx`](results/gsnap_counts.xlsx) on separate tabs.

<!--![](results/assets/screenshot_gsnap_counts.png)-->

(6) Run `multiqc` to see a breakdown of reads by percentage of assigned, mapped/unmapped, etc.

```
multiqc .
mv multiqc_report.html gsnap_multiqc_report.html
```

[gsnap_multiqc_report.html](results/gsnap_multiqc_report.html)

![](results/assets/screenshot_multiqc_gsnap.png)

(7) Added a metadata tab

<!--![](results/assets/screenshot_metadata.png)-->
