# I used https://www.nature.com/articles/nprot.2016.095 (Original nature Protocol paper)
# and the blog entry by Dave Tang (https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/) 
library(genefilter)
library(dplyr)
library(devtools)
BiocManager::install("ballgown")
library(ballgown)
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)
getwd()
# "/Users/Lab/Amy_Toth_discussions/Maize/04_Ballgown/Diff_Expr"
sample_data <- read.csv("~/Amy_Toth_discussions/Maize/S2C.csv")
sample_data
## this is the mapping of the Tissue to sample also called phenotype data
# Sample      Tissue
# 1  SRR1573504    Blade_P7
# 2  SRR1573505    Blade_P7
# 3  SRR1573506    Blade_P7
# 4  SRR1573507   Ligule_P7
# 5  SRR1573508   Ligule_P7
# 6  SRR1573509   Ligule_P7
# 7  SRR1573510   Sheath_P7
# 8  SRR1573511   Sheath_P7
# 9  SRR1573512   Sheath_P7
# 10 SRR1573513       wt_P6
# 11 SRR1573514       wt_P6
# 12 SRR1573515       wt_P6
# 13 SRR1573516      lg1_P6
# 14 SRR1573517      lg1_P6
# 15 SRR1573518      lg1_P6
# 16 SRR1573519  Blade_P7L1
# 17 SRR1573520  Blade_P7L1
# 18 SRR1573521  Blade_P7L1
# 19 SRR1573522 Ligule_P7L1
# 20 SRR1573523 Ligule_P7L1
# 21 SRR1573524 Ligule_P7L1
# 22 SRR1573525 Sheath_P7L1
# 23 SRR1573526 Sheath_P7L1
# 24 SRR1573527 Sheath_P7L1

bg_maize <- ballgown(dataDir = "../../04_Ballgown/",
                    samplePattern = "SRR",
                    pData = sample_data)
class(bg_maize)
# # [1] "ballgown"
# attr(,"package")
# # [1] "ballgown"

# Available methods in class ballgown:

methods(class="ballgown")
# [1] dirs            eexpr           expr            expr<-          geneIDs         geneNames       gexpr           iexpr           indexes        
# [10] indexes<-       mergedDate      pData           pData<-         sampleNames     seqnames        show            structure       subset         
# [19] texpr           transcriptIDs   transcriptNames

dirs(bg_maize) # get all linked directories
head(gexpr(bg_maize),2) 

# FPKM.SRR1573504 FPKM.SRR1573505 FPKM.SRR1573506 FPKM.SRR1573507 FPKM.SRR1573508 FPKM.SRR1573509 FPKM.SRR1573510 FPKM.SRR1573511
# ENSRNA049437471               0               0               0               0               0               0               0               0
# ENSRNA049437473               0               0               0               0               0               0               0               0
# FPKM.SRR1573512 FPKM.SRR1573513 FPKM.SRR1573514 FPKM.SRR1573515 FPKM.SRR1573516 FPKM.SRR1573517 FPKM.SRR1573518 FPKM.SRR1573519
# ENSRNA049437471               0               0               0               0               0               0               0               0
# ENSRNA049437473               0               0               0               0               0               0               0               0
# FPKM.SRR1573520 FPKM.SRR1573521 FPKM.SRR1573522 FPKM.SRR1573523 FPKM.SRR1573524 FPKM.SRR1573525 FPKM.SRR1573526 FPKM.SRR1573527
# ENSRNA049437471               0               0               0               0               0               0               0               0
# ENSRNA049437473               0               0               0               0               0               0               0               0

bg_maize_filtered<-subset(bg_maize,"rowVars(texpr(bg_maize))>1", genomesubset=T)
search()
bg_maize_filtered
# ballgown instance with 47473 transcripts and 24 samples; almost 136K transcripts were found to have low variance and so filtered out
# 
DE_transcripts <- stattest(bg_maize_filtered,
                                feature="transcript",
                                covariate="Tissue",
                                getFC=TRUE, meas="FPKM")
head(DE_transcripts)
# feature id         pval         qval
# 1 transcript  2 2.856799e-06 0.0002857999
# 2 transcript  3 3.795340e-03 0.0573444155
# 3 transcript  4 9.989764e-01 0.9992289906
# 4 transcript  5 9.097779e-02 0.3979534336
# 5 transcript  6 9.790447e-01 0.9863373372
# 6 transcript  7 5.473377e-01 0.8032907650

## Fold change is reported only if we have two sample comparison; so subset the ballgown object to only two Tissue types, viz., wt_P6 and lg1_P6
bg_maize_filtered_P6<-subset(bg_maize_filtered,"Tissue=='wt_P6' | Tissue=='lg1_P6'",genomesubset=F)
DE_transcripts_P6 <- stattest(bg_maize_filtered_P6,
                             feature="transcript",
                             covariate="Tissue",
                             getFC=TRUE, meas="FPKM")
head(DE_transcripts_P6)
# feature id        fc      pval      qval
# 1 transcript  2 0.8941310 0.7753288 0.9601116
# 2 transcript  3 1.0287728 0.9443002 0.9899222
# 3 transcript  4 0.7614619 0.5206759 0.9479091
# 4 transcript  5 0.8027877 0.6936129 0.9479091
# 5 transcript  6 1.4375284 0.8227434 0.9684774
# 6 transcript  7 0.9941932 0.9019003 0.9838098