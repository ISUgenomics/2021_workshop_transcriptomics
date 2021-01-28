# README

## Table of Contents

| File Name | Description |
| -- | -- |
|00_rawdata.md | Getting the raw data |

## Timeline

Introduction of different technologies and their papers in chronological order. I'm probably missing a paper / technology somewhere, so feel free to send suggestions.

**1975** - Sanger Sequencing

* Sanger, F. and Coulson, A.R., 1975. [A rapid method for determining sequences in DNA by primed synthesis with DNA polymerase](https://www.semanticscholar.org/paper/A-rapid-method-for-determining-sequences-in-DNA-by-Sanger-Coulson/8ae9a57dbb5c0a7bbaa13f0ef3fde7a5cd027cd1). Journal of molecular biology, 94(3), pp.441-448.

**1996** - Usually cited as the first mention of nanopore for DNA and RNA seq

* Kasianowicz, J.J., Brandin, E., Branton, D. and Deamer, D.W., 1996. [Characterization of individual polynucleotide molecules using a membrane channel](https://www.semanticscholar.org/paper/Characterization-of-individual-polynucleotide-using-Kasianowicz-Brandin/4f05f0170b5ca6caabd351e2f6f2b1ec10d8a4a8). Proceedings of the National Academy of Sciences, 93(24), pp.13770-13773.

**2008** - First RNAseq paper according to the [Hitzeman et al, 2013](https://www.semanticscholar.org/paper/Genes%2C-behavior-and-next-generation-RNA-sequencing.-Hitzemann-Bottomly/9c05cc4538a272adc0006a2d3e1967b34ce6d1c8) review paper

* Nagalakshmi, U., Wang, Z., Waern, K., Shou, C., Raha, D., Gerstein, M. and Snyder, M., 2008. [The transcriptional landscape of the yeast genome defined by RNA sequencing.](https://www.semanticscholar.org/paper/The-Transcriptional-Landscape-of-the-Yeast-Genome-Nagalakshmi-Wang/447bf5edc72aa9b8c96f841ff11e8f47e89d2ec4) Science, 320(5881), pp.1344-1349.

**2008** - First Brain RNAseq paper according to the [Hitzeman et al, 2013](https://www.semanticscholar.org/paper/Genes%2C-behavior-and-next-generation-RNA-sequencing.-Hitzemann-Bottomly/9c05cc4538a272adc0006a2d3e1967b34ce6d1c8) review paper

* Mortazavi, A., Williams, B.A., McCue, K., Schaeffer, L. and Wold, B., 2008. [Mapping and quantifying mammalian transcriptomes by RNA-Seq](https://www.semanticscholar.org/paper/Mapping-and-quantifying-mammalian-transcriptomes-by-Mortazavi-Williams/ef117c95b92b68b751143155022a5c1a600afe5c). Nature methods, 5(7), pp.621-628.

**2008** - Introducing Illumina Sequencing

* Bentley, D.R., Balasubramanian, S., Swerdlow, H.P., Smith, G.P., Milton, J., Brown, C.G., Hall, K.P., Evers, D.J., Barnes, C.L., Bignell, H.R. and Boutell, J.M., 2008. [Accurate whole human genome sequencing using reversible terminator chemistry](https://www.semanticscholar.org/paper/Accurate-Whole-Human-Genome-Sequencing-using-Bentley-Balasubramanian/6dad16a6941b204c2f5f95d9cda6d0124d5a1a7b). nature, 456(7218), pp.53-59.

**2009** - Introducing PacBio single molecule real-time (SMRT)

* Eid, J., Fehr, A., Gray, J., Luong, K., Lyle, J., Otto, G., Peluso, P., Rank, D., Baybayan, P., Bettman, B. and Bibillo, A., 2009. [Real-time DNA sequencing from single polymerase molecules](https://www.semanticscholar.org/paper/Real-Time-DNA-Sequencing-from-Single-Polymerase-Eid-Fehr/b0588d6e0753e7e3f82c4aa4f1609764179b12de). Science, 323(5910), pp.133-138.

**2010** - EdgeR, usually cited with their 2012 paper

* Robinson, M.D., McCarthy, D.J. and Smyth, G.K., 2010. [edgeR: a Bioconductor package for differential expression analysis of digital gene expression data](https://www.semanticscholar.org/paper/edgeR%3A-a-Bioconductor-package-for-differential-of-Robinson-McCarthy/ec3d71a2fdd01968a6dc638ee261715a0f118c1e). Bioinformatics, 26(1), pp.139-140.
* McCarthy, D.J., Chen, Y. and Smyth, G.K., 2012. [Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation](https://www.semanticscholar.org/paper/Differential-expression-analysis-of-multifactor-to-McCarthy-Chen/571c3ea8cabd16ab0ae7a1a3495d3f3aca918e23). Nucleic acids research, 40(10), pp.4288-4297.

**2012** - QuasiSeq

* Lund, S.P., Nettleton, D., McCarthy, D.J. and Smyth, G.K., 2012. [Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates](https://www.semanticscholar.org/paper/Detecting-Differential-Expression-in-RNA-sequence-Lund-Nettleton/701cae7fb41417989c384bcd8268a620c4669ca4). Statistical applications in genetics and molecular biology, 11(5).

**2013** - STAR

* Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. [STAR: ultrafast universal RNA-seq aligner](https://www.semanticscholar.org/paper/STAR%3A-ultrafast-universal-RNA-seq-aligner-Dobin-Davis/78ce0c149860363bbbe34306c75a4454ad23828d). Bioinformatics, 29(1), pp.15-21.

**2014** - DESeq2

* Love, M.I., Huber, W. and Anders, S., 2014. [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.](https://www.semanticscholar.org/paper/Moderated-estimation-of-fold-change-and-dispersion-Love-Huber/772f5fca88de0f6f38116d73cc32e23efe780a10) Genome biology, 15(12), pp.1-21.

  Hmm...`Love, M.I.` is also a coauthor on a recent RNASeq Review article:

  * Van den Berge, K., Hembach, K.M., Soneson, C., Tiberi, S., Clement, L., Love, M.I., Patro, R. and Robinson, M.D., 2019. [RNA sequencing data: hitchhiker's guide to expression analysis.](https://www.semanticscholar.org/paper/RNA-Sequencing-Data%3A-Hitchhiker's-Guide-to-Analysis-Berge-Hembach/54ffba1e7abd6305f0bd7fcd67d45c202330b25b)

**2014** - featureCounts

* Liao, Y., Smyth, G.K. and Shi, W., 2014. [featureCounts: an efficient general purpose program for assigning sequence reads to genomic features](https://www.semanticscholar.org/paper/featureCounts%3A-an-efficient-general-purpose-program-Liao-Smyth/cdbe8a265ef4bd8350722d2209fc6cc6290da1b3). Bioinformatics, 30(7), pp.923-930.

**2015** - HiSat... which was updated to HiSat2 in 2019?

* Kim, D., Langmead, B. and Salzberg, S.L., 2015. [HISAT: a fast spliced aligner with low memory requirements](https://www.semanticscholar.org/paper/HISAT%3A-a-fast-spliced-aligner-with-low-memory-Kim-Langmead/176490c7bbb20b5a64aca49d9dc3b75bdfd76d67). Nature methods, 12(4), pp.357-360.
* Kim, D., Paggi, J.M., Park, C., Bennett, C. and Salzberg, S.L., 2019. [Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://www.semanticscholar.org/paper/Graph-based-genome-alignment-and-genotyping-with-Kim-Paggi/189b79cda928d58f695cf8323b9ce2196fc22409). Nature biotechnology, 37(8), pp.907-915.

**2018** - Ah...Hitzemann has a recent review paper:

* Hitzemann, R., Darakjian, P., Oberbeck, D., Zheng, C., Walter, N., Iancu, O.D., Searles, R., Harrington, C.A. and McWeeney, S., 2018. [Genes, Behavior, and Next-Generation Sequencing: The First 10 Years](https://www.semanticscholar.org/paper/Genes%2C-behavior%2C-and-next-generation-sequencing%3A-10-Hitzemann-Darakjian/2da210112b07a00d1eb7bece39536e60cb9f8260). In Molecular-Genetic and Statistical Techniques for Behavioral and Neural Research (pp. 289-308). Academic Press.
