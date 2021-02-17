# Install and run rseqc to determine stranding information for maize paired end.

```
#/work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align

ml miniconda3
conda create -n rseqc

source activate rseqc

#this did not finish, using pip
conda install -c bioconda rseq

pip install --user RSeQC

```
