# Literature
This contains literature not included in agenda  that may be helpful to others

1. (Discussed on 24Feb2021) Good overview of basic graph theory: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7004966/
doi: 10.3389/fbioe.2020.00034
  * Shared by Jennifer. Nice organization of paper: Define a network (set of nodes and edges), define types of networks, define centralities, etc

2. MultiQC to assess read alignment (can try this on gsnap output): http://www.bea.ki.se/documents/Intro2RNAseq.pdf

3. FeatureCounts: http://www.bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
  * Understanding featurecounts (with nice diagram of ): https://www.mathworks.com/help/bioinfo/ref/featurecount.html

4. Samtools tutorial: http://quinlanlab.org/tutorials/samtools/samtools.html

5. DESeq2 tutorial: http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink

6. Comparison of RNAseq pipelines:
  * https://www.nature.com/articles/s41598-020-76881-x
  * Consider number of samples, normalization method, batch effect correction, and measure of correlation of expression (Spearman v Pearson) on gene co-expression networks: https://www.biorxiv.org/content/10.1101/2021.03.11.435043v1
    * Pearson: based on raw values (sensitive to extreme values), better for medium or large datasets (> 50)
    * Spearman: based on ranked values, better for small datasets (< 30)
    * ComBat to treat batch effects (https://bioconductor.org/packages/release/bioc/html/sva.html)
    * UQ (upper quartile) normalization
