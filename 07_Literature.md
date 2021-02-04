# Literature related to RNA-Seq analysis




| Paper (and link to the paper) | year published | dataset | Alignment methods used | DGE identifying methods used | Any alignment tool suggested to work better in this paper | Any DGE identifier tool suggested to perform better in this paper | Is it supported by experimental data ? | person found the paper | Notes |
| --- | ---| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| [Evaluation of Seven Different RNA-Seq AlignmentTools Based on Experimental Data from the ModelPlantArabidopsis thaliana](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjS8Nnw7s7uAhX7Ap0JHbwrB3gQFjABegQIBBAC&url=https%3A%2F%2Fwww.mdpi.com%2F1422-0067%2F21%2F5%2F1720%2Fpdf&usg=AOvVaw0wKcqBj60sas8h5syWD39b)| 2020 | Arabidopsis thaliana | bwa,  CLC Genomics Workbench,  HISAT2,   kallisto,  RSEM, salmon and  STAR | DESeq2, CLC  |STAR is better in alignment (higher tolerance for soft-clipping) but the final DGE results are similar between methods|CLC produced 50% more DE genes compare to DESeq2, can't confirm if this is better |  yes but only for alignment |   bwa, salmon and kallisto, using the transcriptomic reference, identified less genes. This difference is due to the presence of non-coding RNAs such as transfer RNAs (tRNA) and micro RNAs (miRNA) in the genomic reference, which are absent from the transcriptomic reference -  transcripts with less than five counts were filtered  | Maryam |  
