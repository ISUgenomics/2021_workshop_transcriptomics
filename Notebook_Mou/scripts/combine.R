#! /usr/bin/env Rscript
# Auth: Jennifer Chang & Mou
# Date: 2021/03/04
# Desc: Combine featureCounts output (1st and last column) files for Bee and Maize. The text files (*.genecounts.txt) generated from this script will be used for DESeq2.

# === Load Libraries
library(tidyverse)
library(magrittr)
library(readxl)

###### BEE #######

# === Get list of featureCount output files
dir_org="~/Desktop/bee/"        # counts are in a "bee" or "maize" subdirectory
featureCount_files <- list.files(path = dir_org, pattern = "*genecounts.txt$", full.names = TRUE)

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" ) %>%
  select(Geneid, ends_with(".bam")) %>%              # Get 1st and last column (column was named after bam file)
  pivot_longer(cols=ends_with(".bam")) %>%           # Melt data (tidy data)
  mutate(
    name = gsub(".Aligned.sortedByCoord.out.bam", "", name)        # No longer need the bam extension, easier to read
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(Geneid, ends_with(".bam")) %>%
    pivot_longer(cols=ends_with(".bam")) %>%
    mutate(
      name = gsub(".Aligned.sortedByCoord.out.bam", "", name)
    )
  data = rbind(data, temp)
}

# === Convert to excel like data (wider)
wide_data <- data %>%
  pivot_wider(id_cols=Geneid)

# === Save tab delimited file (smaller file size)
write_delim(wide_data,
            paste(dir_org, "bee.genecounts.txt"),
            delim="\t")

# === Save Excel file (can be easier to work with)
writexl::write_xlsx(wide_data,
                    path=paste(dir_org, "bee.genecounts.xlsx"))


###### MAIZE #######

# === Get list of featureCount output files
dir_org="~/Desktop/maize/"        # counts are in a "bee" or "maize" subdirectory
featureCount_files <- list.files(path = dir_org, pattern = "*genecounts.txt$", full.names = TRUE)

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" ) %>%
  select(Geneid, ends_with(".bam")) %>%              # Get 1st and last column (column was named after bam file)
  pivot_longer(cols=ends_with(".bam")) %>%           # Melt data (tidy data)
  mutate(
    name = gsub(".aligned.out.bam", "", name)        # No longer need the bam extension, easier to read
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(Geneid, ends_with(".bam")) %>%
    pivot_longer(cols=ends_with(".bam")) %>%
    mutate(
      name = gsub(".aligned.out.bam", "", name)
    )
  data = rbind(data, temp)
}

# === Convert to excel like data (wider)
wide_data <- data %>%
  pivot_wider(id_cols=Geneid)

# === Save tab delimited file (smaller file size)
write_delim(wide_data,
            paste(dir_org, "maize.genecounts.txt"),
            delim="\t")

# === Save Excel file (can be easier to work with)
writexl::write_xlsx(wide_data,
                    path=paste(dir_org, "maize.genecounts.xlsx"))
