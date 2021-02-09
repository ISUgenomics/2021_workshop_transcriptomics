# 2021-01-25 Day 1

## [Genes, behavior and next‐generation RNA sequencing](https://onlinelibrary.wiley.com/doi/full/10.1111/gbb.12007) 2013
 	Hitzemann, R., Bottomly, D., Darakjian, P., Walter, N., Iancu, O., Searles, R., Wilmot, B. and McWeeney, S., 2013. Genes, behavior and next‐generation RNA sequencing. Genes, Brain and Behavior, 12(1), pp.1-12.


## [RNA sequencing and analysis.](http://cshprotocols.cshlp.org/content/2015/11/pdb.top084970.full) 2015
Kukurba, K.R. and Montgomery, S.B., 2015.  Cold Spring Harbor Protocols, 2015(11), pp.pdb-top084970.



# Nanopore RNA-Seq

### https://www.biorxiv.org/content/10.1101/2021.01.22.427687v1

### [A comprehensive examination of Nanopore native RNA sequencing for characterization of complex transcriptomes](https://www.nature.com/articles/s41467-019-11272-z)
A platform for highly parallel direct sequencing of native RNA strands was recently described by Oxford Nanopore Technologies, but despite initial efforts it remains crucial to further investigate the technology for quantification of complex transcriptomes. Here we undertake native RNA sequencing of polyA + RNA from two human cell lines, analysing ~5.2 million aligned native RNA reads. To enable informative comparisons, we also perform relevant ONT direct cDNA- and Illumina-sequencing. We find that while native RNA sequencing does enable some of the anticipated advantages, key unexpected aspects currently hamper its performance, most notably the quite frequent inability to obtain full-length transcripts from single reads, as well as difficulties to unambiguously infer their true transcript of origin. While characterising issues that need to be addressed when investigating more complex transcriptomes, our study highlights that with some defined improvements, native RNA sequencing could be an important addition to the mammalian transcriptomics toolbox.


### [Transcriptome variation in human tissues revealed by long-read sequencing](https://www.biorxiv.org/content/10.1101/2021.01.22.427687v1)

Regulation of transcript structure generates transcript diversity and plays an important role in human disease. The advent of long-read sequencing technologies offers the opportunity to study the role of genetic variation in transcript structure. In this paper, we present a large human long-read RNA-seq dataset using the Oxford Nanopore Technologies platform from 88 samples from GTEx tissues and cell lines, complementing the GTEx resource. We identified just under 100,000 new transcripts for annotated genes, and validated the protein expression of a similar proportion of novel and annotated transcripts. We developed a new computational package, LORALS, to analyze genetic effects of rare and common variants on the transcriptome via allele-specific analysis of long reads. We called allele-specific expression and transcript structure events, providing novel insights into the specific transcript alterations caused by common and rare genetic variants and highlighting the resolution gained from long-read data. We were able to perturb transcript structure upon knockdown of PTBP1, an RNA binding protein that mediates splicing, thereby finding genetic regulatory effects that are modified by the cellular environment. Finally, we use this dataset to enhance variant interpretation and study rare variants leading to aberrant splicing patterns.


# Feb 3, 2021 Week 2 Day 1

### [Aligning the Aligners: Comparison of RNA Sequencing Data Alignment and Gene Expression Quantification Tools for Clinical Breast Cancer Research](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6617288/)

Really interesting read about a polyA tail bias that HiSat2 has over STAR in alignment that could lead to false positive DGE genes.

### [Compare and Contrast of Differential Gene Expression Software Packages of RNA-Seq](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8947662)

This paper compared the different ways to perform differential Gene expression analysis: DESeq2, DEGSeq, edgeR and cuffDiff2.  It tries to make the case for cuffDiff2 but I am not convinced as their criteria for being better appears to be a lower DGE count in the list.  We know from the previous paper that these inflated counts may be due to polyA tail in alignment with bowtie2/hisat2 which was used as the aligner in this paper.

### [RNA-Seq differential expression analysis: An extended review and a software tool](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190152)
