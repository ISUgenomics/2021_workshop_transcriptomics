# 06 Kallisto

(1) Index reference transcriptome

```
module load kallisto

GENOME_CDNA=02_Genome/GCF_000188095.3_BIMP_2.2_rna.fna.gz
G_NAME=$(basename $GENOME_CDNA) |sed 's/.fna.gz//g'

kallisto index -i ${G_NAME}.idx ${GENOME_CDNA}
```

(2) Quantify gene counts

```
PROC=16

for READ in 00_Raw-Data/*1.fastq.gz
do
  READNAME=$(basename ${FILE}) |sed 's/.fastq.gz//g'
  READ2=$(echo $READ | sed 's/1.fastq.gz/2.fastq.gz/g')
  
  kallisto quant -i ${G_NAME}.idx -o ${READNAME}_quant -b 20 -t ${PROC} ${READ} ${READ2}
done
```

Counts are in `${READNAME}_quant/abundance.tsv`

```
target_id       length  eff_length      est_counts      tpm
NM_001111367.2  1136    983.207 10149   174.858
NM_001111369.2  1884    1731.21 5219.47 51.0721
NM_001111370.3  2299    2146.21 0       0
NM_001111373.2  1660    1507.21 7208.67 81.0195
NM_001111374.2  791     638.207 0       0
...
```

(3) Combine the Kallisto counts (`est_counts` and `tsv`) into one file, by modifying the R script from GSNAP

```
#! /usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(readxl)

# === Get list of featureCount output files
(featureCount_files <- list.files(path = "03_Kallisto", pattern = "*_quant$", full.names = TRUE) %>% 
    paste(., "/abundance.tsv", sep=""))

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" ) %>%
  select(target_id, est_counts, tpm) %>%             # est_counts or tpm
  pivot_longer(cols=c(est_counts, tpm)) %>% 
  mutate(
    sample=featureCount_files[1] %>% gsub("03_Kallisto/","",.) %>% gsub("_quant/abundance.tsv","", .)
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(target_id, est_counts, tpm) %>%             # est_counts or tpm
    pivot_longer(cols=c(est_counts, tpm)) %>% 
    mutate(
      sample=count_file %>% gsub("03_Kallisto/","",.) %>% gsub("_quant/abundance.tsv","", .)
    )
  data = rbind(data, temp)
}

# === Convert to excel like data (wider)
data_est = subset(data, name=="est_counts")
data_tpm = subset(data, name=="tpm")

wide_est = data_est %>%
  pivot_wider(id_cols=target_id, names_from=sample, values_from=value )

wide_tpm = data_tpm %>%
  pivot_wider(id_cols=target_id, names_from=sample, values_from=value )
  

# === Save tab delimited file (smaller file size)
write_delim(wide_est, "kallisto_est_counts.txt", delim="\t")
write_delim(wide_tpm, "kallisto_tpm.txt", delim="\t")
```

Then combining `kallisto_est_counts.txt` and `kallisto_tpm.txt` in one `kallisto_counts.xlsx` with separate tabs.

<!--![](results/assets/screenshot_kallisto_counts.png)-->

View [kallisto_counts.xlsx](results/kallisto_counts.xlsx).





