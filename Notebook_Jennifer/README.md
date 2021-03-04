# README

## Table of Contents

| File Name | Description | Year of Development |
|:--|:--|:-:|
|[00_data.md](00_data.md) | Getting the raw data onto Atlas HPC| |
|[01\_align\_gsnap.md](01_align_gsnap.md) | GSNAP alignment| 2010 |
|[01\_align\_salmon.md](01_align_salmon.md) | Salmon alignment... in progress, ready to test | 2017 |
|[01\_align\_hisat2.md](01_align_hisat2.md) | HISAT2 alignment (with ballgown)... in progress| 2015 |
|[02\_DE\_analysis.md](02_DE_analysis.md) | Installing DESeq2 in R ... in progress | 2014 |
|[03\_Normalization.md](02_Normalization.md) | placeholder for normalization methods |  |
|[04\_WGCNA.md](04_WGCNA.md) | WGCNA (networks) maize | 2008 |
|[04\_WGCNA_bee.md](04_WGCNA_bee.md) | WGCNA (networks) bee | 2008 |

## Timeline

Introduction of different technologies and their academic papers in chronological order.

**1975** - Sanger Sequencing - type:seqtech

* Sanger, F. and Coulson, A.R., 1975. [A rapid method for determining sequences in DNA by primed synthesis with DNA polymerase](https://api.semanticscholar.org/CorpusID:30906971). Journal of molecular biology, 94(3), pp.441-448.

**1991** - Introducing microarrays - type:expressiontech

* Fodor, S.P., Read, J.L., Pirrung, M.C., Stryer, L., Lu, A.T. and Solas, D., 1991. [Light-directed, spatially addressable parallel chemical synthesis](https://api.semanticscholar.org/CorpusID:14754629). science, 251(4995), pp.767-773.

**1996** - Usually cited as the first mention of nanopore for DNA and RNA seq - type:seqtech

* Kasianowicz, J.J., Brandin, E., Branton, D. and Deamer, D.W., 1996. [Characterization of individual polynucleotide molecules using a membrane channel](https://api.semanticscholar.org/CorpusID:285178). Proceedings of the National Academy of Sciences, 93(24), pp.13770-13773.

**1999** - KEGG - type: biological pathway networks

* Ogata, H., Goto, S., Sato, K., Fujibuchi, W., Bono, H. and Kanehisa, M., 1999. [KEGG: Kyoto encyclopedia of genes and genomes](https://api.semanticscholar.org/CorpusID:7449269). Nucleic acids research, 27(1), pp.29-34.

**2000** - Gene Ontology - type: semantic networks

* Ashburner, M., Ball, C.A., Blake, J.A., Botstein, D., Butler, H., Cherry, J.M., Davis, A.P., Dolinski, K., Dwight, S.S., Eppig, J.T. and Harris, M.A., 2000. [Gene ontology: tool for the unification of biology](https://api.semanticscholar.org/CorpusID:10718909). Nature genetics, 25(1), pp.25-29.

**2001** - GraphViz - type:networks

* Ellson, J., Gansner, E., Koutsofios, L., North, S.C. and Woodhull, G., 2001, September. [Graphviz—open source graph drawing tools](https://api.semanticscholar.org/CorpusID:35525078). In International Symposium on Graph Drawing (pp. 483-484). Springer, Berlin, Heidelberg.

**2003** - Cytoscape - type:networks

* Shannon, P., Markiel, A., Ozier, O., Baliga, N.S., Wang, J.T., Ramage, D., Amin, N., Schwikowski, B. and Ideker, T., 2003. [Cytoscape: a software environment for integrated models of biomolecular interaction networks](https://api.semanticscholar.org/CorpusID:15588516). Genome research, 13(11), pp.2498-2504.

**2008** - First RNAseq paper according to the [Hitzeman et al, 2013](https://api.semanticscholar.org/CorpusID:32346093) review paper

* Nagalakshmi, U., Wang, Z., Waern, K., Shou, C., Raha, D., Gerstein, M. and Snyder, M., 2008. [The transcriptional landscape of the yeast genome defined by RNA sequencing.](https://api.semanticscholar.org/CorpusID:206513052) Science, 320(5881), pp.1344-1349.

**2008** - First Brain RNAseq paper according to the [Hitzeman et al, 2013](https://api.semanticscholar.org/CorpusID:32346093) review paper

* Mortazavi, A., Williams, B.A., McCue, K., Schaeffer, L. and Wold, B., 2008. [Mapping and quantifying mammalian transcriptomes by RNA-Seq](https://api.semanticscholar.org/CorpusID:205418589). Nature methods, 5(7), pp.621-628.

**2008** - WGCNA, theory paper from 2005 - type:networks

* Langfelder, P. and Horvath, S., 2008. [WGCNA: an R package for weighted correlation network analysis](https://api.semanticscholar.org/CorpusID:206970636). BMC bioinformatics, 9(1), pp.1-13.
* Zhang, B. and Horvath, S., 2005. [A general framework for weighted gene co-expression network analysis](https://api.semanticscholar.org/CorpusID:7756201). Statistical applications in genetics and molecular biology, 4(1).

**2008** - Arena3D, with updates in 2012, and 2020 (web interface?) - type:networks

* Pavlopoulos, G.A., O'Donoghue, S.I., Satagopam, V.P., Soldatos, T.G., Pafilis, E. and Schneider, R., 2008. [Arena3D: visualization of biological networks in 3D](https://api.semanticscholar.org/CorpusID:9118716). BMC systems biology, 2(1), pp.1-7.
* Secrier, M., Pavlopoulos, G.A., Aerts, J. and Schneider, R., 2012. [Arena3D: visualizing time-driven phenotypic differences in biological systems](https://api.semanticscholar.org/CorpusID:5068898). BMC bioinformatics, 13(1), pp.1-11.
* Karatzas, E., Baltoumas, F.A., Panayiotou, N.A., Schneider, R. and Pavlopoulos, G.A., 2020. [Arena3Dweb: Interactive 3D visualization of multilayered networks](https://api.semanticscholar.org/CorpusID:227172915). bioRxiv.

**2008** - Introducing Illumina Sequencing - type:seqtech

* Bentley, D.R., Balasubramanian, S., Swerdlow, H.P., Smith, G.P., Milton, J., Brown, C.G., Hall, K.P., Evers, D.J., Barnes, C.L., Bignell, H.R. and Boutell, J.M., 2008. [Accurate whole human genome sequencing using reversible terminator chemistry](https://api.semanticscholar.org/CorpusID:4417841). nature, 456(7218), pp.53-59.

**2009** - Introducing PacBio single molecule real-time (SMRT) - type:seqtech

* Eid, J., Fehr, A., Gray, J., Luong, K., Lyle, J., Otto, G., Peluso, P., Rank, D., Baybayan, P., Bettman, B. and Bibillo, A., 2009. [Real-time DNA sequencing from single polymerase molecules](https://api.semanticscholar.org/CorpusID:54488479). Science, 323(5910), pp.133-138.

**2009** - Bowtie, updated to Bowtie2 in 2012 - type:aligner

* Langmead, B., Trapnell, C., Pop, M. and Salzberg, S.L., 2009. [Ultrafast and memory-efficient alignment of short DNA sequences to the human genome](https://api.semanticscholar.org/CorpusID:5057). Genome biology, 10(3), pp.1-10.
* Langmead, B. and Salzberg, S.L., 2012. [Fast gapped-read alignment with Bowtie 2](https://api.semanticscholar.org/CorpusID:205420407). Nature methods, 9(4), p.357.

* Bowtie2 GitHub - [https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)

**2009** - TopHat, updated to TopHat2 in 2013 - type:aligner

* Trapnell, C., Pachter, L. and Salzberg, S.L., 2009. [TopHat: discovering splice junctions with RNA-Seq](https://api.semanticscholar.org/CorpusID:1866358). Bioinformatics, 25(9), pp.1105-1111.
* Kim, D., Pertea, G., Trapnell, C., Pimentel, H., Kelley, R. and Salzberg, S.L., 2013. [TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions](https://api.semanticscholar.org/CorpusID:1663277). Genome biology, 14(4), pp.1-13.

* TopHat Website - [http://ccb.jhu.edu/software/tophat/index.shtml](http://ccb.jhu.edu/software/tophat/index.shtml)

**2009** - BWA - type:aligner

* Li, H. and Durbin, R., 2009. [Fast and accurate short read alignment with Burrows–Wheeler transform](https://api.semanticscholar.org/CorpusID:3662132). bioinformatics, 25(14), pp.1754-1760.

**2010** - EdgeR, usually cited with their 2012 paper - type:DEanalysis

* Robinson, M.D., McCarthy, D.J. and Smyth, G.K., 2010. [edgeR: a Bioconductor package for differential expression analysis of digital gene expression data](https://api.semanticscholar.org/CorpusID:1481014). Bioinformatics, 26(1), pp.139-140.
* McCarthy, D.J., Chen, Y. and Smyth, G.K., 2012. [Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation](https://api.semanticscholar.org/CorpusID:14052240). Nucleic acids research, 40(10), pp.4288-4297.

* EdgeR Bioconductor - [https://bioconductor.org/packages/release/bioc/html/edgeR.html](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

**2010** - baySeq - type:DEanalysis

* Hardcastle, T.J. and Kelly, K.A., 2010. [baySeq: empirical Bayesian methods for identifying differential expression in sequence count data](https://api.semanticscholar.org/CorpusID:1133459). BMC bioinformatics, 11(1), pp.1-14.

* baySeq Bioconductor - [https://bioconductor.org/packages/release/bioc/html/baySeq.html](https://bioconductor.org/packages/release/bioc/html/baySeq.html)

**2010** - GSNAP - type:aligner

* Wu, T.D. and Nacu, S., 2010. [Fast and SNP-tolerant detection of complex variants and splicing in short reads](https://api.semanticscholar.org/CorpusID:15689019). Bioinformatics, 26(7), pp.873-881.

**2012** - QuasiSeq - type:DEanalysis

* Lund, S.P., Nettleton, D., McCarthy, D.J. and Smyth, G.K., 2012. [Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates](https://api.semanticscholar.org/CorpusID:8580642). Statistical applications in genetics and molecular biology, 11(5).

* QuasiSeq CRAN - [https://cran.r-project.org/web/packages/QuasiSeq/index.html](https://cran.r-project.org/web/packages/QuasiSeq/index.html)

**2013** - STAR - type:aligner

* Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. [STAR: ultrafast universal RNA-seq aligner](https://api.semanticscholar.org/CorpusID:11244195). Bioinformatics, 29(1), pp.15-21.

* STAR GitHub - [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)

<!-- This paper goes with genome assembly software

**2013** - MaSuRCA genome assembly

* Zimin, A.V., Marçais, G., Puiu, D., Roberts, M., Salzberg, S.L. and Yorke, J.A., 2013. [The MaSuRCA genome assembler](https://www.semanticscholar.org/paper/The-MaSuRCA-genome-assembler-Zimin-Mar%C3%A7ais/476273d8efa79c3067e55bc5b556e34c88225491). Bioinformatics, 29(21), pp.2669-2677.
-->

**2013** - CUFFDIFF2 - type:DEanalysis

* Trapnell, C., Hendrickson, D.G., Sauvageau, M., Goff, L., Rinn, J.L. and Pachter, L., 2013. [Differential analysis of gene regulation at transcript resolution with RNA-seq](https://api.semanticscholar.org/CorpusID:9253369). Nature biotechnology, 31(1), pp.46-53.

* CUFFDIFF website - [http://cole-trapnell-lab.github.io/cufflinks/getting_started/](http://cole-trapnell-lab.github.io/cufflinks/getting_started/)

**2014** - DESeq2 - type:DEanalysis

* Love, M.I., Huber, W. and Anders, S., 2014. [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.](https://api.semanticscholar.org/CorpusID:8620592) Genome biology, 15(12), pp.1-21.

  Hmm...`Love, M.I.` is also a coauthor on a recent RNASeq Review article:

  `Love, M.I.` (DESeq developer) and `Patro,R ` (Salmon developer) co-authored a paper together, and also a review. 

  * Love, M.I., Soneson, C. and Patro, R., 2018. [Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification](https://api.semanticscholar.org/CorpusID:53024683). F1000Research, 7.
  * Van den Berge, K., Hembach, K.M., Soneson, C., Tiberi, S., Clement, L., Love, M.I., Patro, R. and Robinson, M.D., 2019. [RNA sequencing data: hitchhiker's guide to expression analysis.](https://api.semanticscholar.org/CorpusID:226923140)

* DESeq2 Bioconductor - [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

**2014** - featureCounts - type:count?

* Liao, Y., Smyth, G.K. and Shi, W., 2014. [featureCounts: an efficient general purpose program for assigning sequence reads to genomic features](https://api.semanticscholar.org/CorpusID:15960459). Bioinformatics, 30(7), pp.923-930.

* featureCount Website - [http://subread.sourceforge.net/](http://subread.sourceforge.net/)

**2014** - Ballgown, originally on bioconductor in 2014 but the citation always moves to the latest version - type:DEanalysis

* Fu J, Frazee AC, Collado-Torres L, Jaffe AE, Leek JT (2020). [ballgown: Flexible, isoform-level differential expression analysis](https://api.semanticscholar.org/CorpusID:17545826). R package version 2.22.0.

* Ballgown Bioconductor - [https://www.bioconductor.org/packages/release/bioc/html/ballgown.html](https://www.bioconductor.org/packages/release/bioc/html/ballgown.html)

**2015** - NOISeq R - type:DEanalysis

* Tarazona, S., Furió-Tarí, P., Turrà, D., Pietro, A.D., Nueda, M.J., Ferrer, A. and Conesa, A., 2015. [Data quality aware analysis of differential expression in RNA-seq with NOISeq R/Bioc package](https://api.semanticscholar.org/CorpusID:16206150). Nucleic acids research, 43(21), pp.e140-e140.

* NOISeq Bioconductor - [https://bioconductor.org/packages/release/bioc/html/NOISeq.html](https://bioconductor.org/packages/release/bioc/html/NOISeq.html)

**2015** - HiSat... which was updated to HiSat2 in 2019? - type:aligner

* Kim, D., Langmead, B. and Salzberg, S.L., 2015. [HISAT: a fast spliced aligner with low memory requirements](https://api.semanticscholar.org/CorpusID:15738431). Nature methods, 12(4), pp.357-360.
* Kim, D., Paggi, J.M., Park, C., Bennett, C. and Salzberg, S.L., 2019. [Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://api.semanticscholar.org/CorpusID:199388685). Nature biotechnology, 37(8), pp.907-915.

* HiSat Website - [http://daehwankimlab.github.io/hisat2/](http://daehwankimlab.github.io/hisat2/)

**2016** - Mango Graph Studio (altho this feels like self promotion), not restricted to biological networks, applied to literature mining in Cavatica - type:networks

* Chang, J., Cho, H. and Chou, H.H., 2016. [Mango: combining and analyzing heterogeneous biological networks](https://api.semanticscholar.org/CorpusID:15220779). BioData mining, 9(1), pp.1-14.
* Chang, J. and Chou, H.H., 2017, November. [Cavatica: A pipeline for identifying author adoption trends among software or methods](https://api.semanticscholar.org/CorpusID:659763). In 2017 IEEE International Conference on Bioinformatics and Biomedicine (BIBM) (pp. 2145-2150). IEEE.

**2016** - Kallisto - type:workflow

* Bray, N.L., Pimentel, H., Melsted, P. and Pachter, L., 2016. [Near-optimal probabilistic RNA-seq quantification](https://api.semanticscholar.org/CorpusID:205282743). Nature biotechnology, 34(5), pp.525-527.

* Kallisto Website - [https://pachterlab.github.io/kallisto/about](https://pachterlab.github.io/kallisto/about)

**2017** - Salmon, with discussion of methods in 2019/20ish - type:workflow

* Patro, R., Duggal, G., Love, M.I., Irizarry, R.A. and Kingsford, C., 2017. [Salmon provides fast and bias-aware quantification of transcript expression](https://api.semanticscholar.org/CorpusID:195671818). Nature methods, 14(4), pp.417-419.
* Srivastava, A., Malik, L., Sarkar, H., Zakeri, M., Almodaresi, F., Soneson, C., Love, M.I., Kingsford, C. and Patro, R., 2020. [Alignment and mapping methodology influence transcript abundance estimation](https://api.semanticscholar.org/CorpusID:195435380). Genome biology, 21(1), pp.1-29.

  Hmm... seems to be some strange discussion on Salmon vs Kallisto
  
  * [https://liorpachter.wordpress.com/2017/08/02/how-not-to-perform-a-differential-expression-analysis-or-science/](https://liorpachter.wordpress.com/2017/08/02/how-not-to-perform-a-differential-expression-analysis-or-science/)

* Salmon Website - [https://combine-lab.github.io/salmon/](https://combine-lab.github.io/salmon/)


**2018** - Ah...Hitzemann has a recent review paper:

* Hitzemann, R., Darakjian, P., Oberbeck, D., Zheng, C., Walter, N., Iancu, O.D., Searles, R., Harrington, C.A. and McWeeney, S., 2018. [Genes, Behavior, and Next-Generation Sequencing: The First 10 Years](https://api.semanticscholar.org/CorpusID:91255778). In Molecular-Genetic and Statistical Techniques for Behavioral and Neural Research (pp. 289-308). Academic Press.

**2020** - CoCoCoNet - type:network analysis

* Lee, J., Shah, M., Ballouz, S., Crow, M. and Gillis, J., 2020. [CoCoCoNet: conserved and comparative co-expression across a diverse set of species](https://api.semanticscholar.org/CorpusID:218600760). Nucleic acids research, 48(W1), pp.W566-W571.

## Brainstorming questions

### What are the earliest aligner and DE analysis programs?

Maybe Bowtie and TopHat (earliest papers at 2009), and both were updated in 2012.

### What are the most recently developed aligner and DE analysis programs?

Maybe HiSat2 (2019), Salmon (2017), and NOISeq R (2015)? The Salmon and DESeq people seem to be working together.
