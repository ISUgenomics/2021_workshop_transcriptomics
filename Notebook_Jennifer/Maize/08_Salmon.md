# 08 Salmon

(1) Index reference transcriptome

```
module load salmon

GENOME_CDNA=02_Genome/GCF_000188095.3_BIMP_2.2_rna.fna.gz
G_NAME=$(basename $GENOME_CDNA) |sed 's/.fna.gz//g'

salmon index -i ${G_NAME} -t ${GENOME}
```

(2) Quantify gene counts

```
PROC=16

for READ in 00_Raw-Data/*1.fastq.gz
do
  READNAME=$(basename ${FILE}) |sed 's/.fastq.gz//g'
  READ2=$(echo $READ | sed 's/1.fastq.gz/2.fastq.gz/g')

  salmon quant -l A -p ${PROC} --validateMappings -i ${G_NAME} -1 ${READ} -2 ${READ2} -o ${READNAME}_quant
done
```

The gene counts are in ``${READNAME}_quant/quant.sf`

```
Name    Length  EffectiveLength TPM     NumReads
NM_001111367.2  1136    976.352 189.327741      10168.000
NM_001111369.2  1884    1708.352        52.197159       4905.000
NM_001111370.3  2299    2122.352        0.000000        0.000
NM_001111373.2  1660    1500.352        82.370784       6798.000
NM_001111374.2  791     631.365 0.000000        0.000
NM_001111375.2  843     683.362 202.876251      7626.000
...
```

(3) Run `multiqc`

```
cd 03_Salmon
multiqc .
mv multiqc_report.html salmon_multiqc_report.html
```

View [salmon_multiqc_report.html](results/salmon_multiqc_report.html)

(4) Combine counts information into one file, by modifying the script for GSNAP

```
#! /usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(readxl)

# === Get list of featureCount output files
(featureCount_files <- list.files(path = "03_Salmon", pattern = "*_quant$", full.names = TRUE) %>% 
    paste(., "/quant.sf", sep=""))

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" ) %>%
  select(Name, TPM, NumReads) %>%             # est_counts or tpm
  pivot_longer(cols=c(TPM, NumReads)) %>% 
  mutate(
    sample=featureCount_files[1] %>% gsub("03_Salmon/","",.) %>% gsub("_quant/quant.sf","", .)
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(Name, TPM, NumReads) %>%             # est_counts or tpm
    pivot_longer(cols=c(TPM, NumReads)) %>% 
    mutate(
      sample=count_file %>% gsub("03_Salmon/","",.) %>% gsub("_quant/quant.sf","", .)
    )
  data = rbind(data, temp)
}

# === Convert to excel like data (wider)
data_numreads = subset(data, name=="NumReads")
data_tpm = subset(data, name=="TPM")

wide_est = data_numreads %>%
  pivot_wider(id_cols=Name, names_from=sample, values_from=value )

wide_tpm = data_tpm %>%
  pivot_wider(id_cols=Name, names_from=sample, values_from=value )
  

# === Save tab delimited file (smaller file size)
write_delim(wide_est, "salmon_numreads.txt", delim="\t")
write_delim(wide_tpm, "salmon_tpm.txt", delim="\t")
```

Combine `salmon_numreads.txt` and `salmon_tpm.txt` into an excel file, on separate tabs.

<!--![](results/assets/screenshot_salmon_counts.png)-->

View [salmon_counts.xlsx](results/salmon_counts.xlsx)
