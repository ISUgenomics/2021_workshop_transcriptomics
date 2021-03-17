# Results


Comparing different alignments and DGE identifying Methods
|dataset 	|Alignment method |method of counting 	|DEG identifying method 	|Who is working on this problem|# of DE genes|
| --- | --- | --- |---| --- | --- |
|Bee 	|Star | [featureCounts](./data/counts/Maryam-STAR-featureCounts-At_count.txt)|DESeq2 	|Maryam and Ambi| |
|Maize 	|Star | featureCounts|DESeq2 	|Maryam and Ambi| |
|Bee 	|Hisat2 |	|DESeq2 	|Rick and Alex||
|Maize 	|Hisat2 |	|DESeq2 	|Rick and Alex||
|Bee 	|gsnap |[featureCounts](Notebook_Jennifer/bee_genecounts.txt), [featureCounts2](Notebook_Mou/results/bee.genecounts.out.txt)	|DESeq2 	|Jennifer and Kathy||
|Maize 	|gsnap |[featureCounts](Notebook_Jennifer/maize_genecounts.txt), [featureCounts2](Notebook_Mou/results/maize.genecounts.out.txt)	|DESeq2 	|Jennifer and Kathy||
|Bee 	|gsnap |	|EdgeR 	|Severin||
|Maize 	|gsnap |	|EdgeR 	|Severin||
|Bee 	|Star |	|EdgeR 	|Sathesh and Katie| |
|Maize	|Star |	|EdgeR 	|Sathesh and Katie| |
|Bee 	|Hisat2 |[featureCounts](https://github.com/ISUgenomics/2021_workshop_transcriptomics/blob/main/data/counts/Ryan_Bee_count_table.txt)	|EdgeR |	Ryan and Siva||
|Maize	|Hisat2 |	|EdgeR |	Ryan and Siva||
|Bee 	|Hisat2 |	|Stringtie/Ballgown 	|Siva and Jennifer||
|Maize	|Hisat2 |	|Stringtie/Ballgown 	|Siva and Jennifer||
|Bee 	|Kallisto |	|EdgeR 	|Severin/Kyle||
|Maize	|Kallisto |	|EdgeR 	|Severin/Kyle||
|Bee 	|Salmon| 	|EdgeR 	|Jennifer/Kyle||
|Maize 	|Salmon| 	|EdgeR 	|Jennifer/Kyle||
