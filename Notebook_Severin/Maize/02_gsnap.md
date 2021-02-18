# Align reads

* Feb 16, 2021
* /work/gif/TranscriptomicsWorkshop/severin/Maize/02_gsnap/

# Maize

### Softlink the files

```
ls /work/gif/TranscriptomicsWorkshop/Maize/*fastq | more | xargs -I xx ln -s xx
```

### Download the Maize genome

* /work/gif/TranscriptomicsWorkshop/severin/Maize/02_gsnap/assembly
Download and combine chromosomes

```
# V4
#wget ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.*gz
zcat *.gz > Zea_All.fasta
# V5
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
```

### Build the Maize database for gsnap/gmap

*  /work/gif/TranscriptomicsWorkshop/severin/Maize/02_gsnap

```
gmap_build -D Zea -d B73 Zea_All.fasta
gmap_build -D Zea -d NAMV5 GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
```


### Run GSNAP using Nextflow

I am going to use my nextflow parallel workflow to run this.


#### Create a gsnap script

To make it easier to understand I am going to put the command in a simple bash script.


This script does the following
* loads gsnap module
* creates varaibles for the GMAPDB, DBNAME, Reads1 and Reads2 and output file name
* executes gsnap

####  gsnapScript.sh
```
#!/bin/bash

module load gmap-gsnap/2019-05-12-zjqshxf
export GMAPDB=/work/gif/TranscriptomicsWorkshop/severin/Maize/02_gsnap/Zea/NAMV5/
DB_NAME="NAMV5"
FILE1="$1"
FILE2="$2"
OUTFILE=$(basename ${FILE1} | sed 's/_1.fastq$//g')
# Note: "-N" option for detecting novel splice sites, remove if not needed (0=OFF; 1=ON)
gsnap -d ${DB_NAME} -N 1 -t 8 -B 4 -m 5 --input-buffer-size=1000000 --output-buffer-size=1000000 -A sam --split-output=${DB_NAME}_${OUTFILE} ${FILE1} ${FILE2}
```

#### run the script using nextflow parallel workflow

Parallel workflow requires the input files command which will be executed and placed in a queue that will be fed to the script parameter which in our case is the script above.

```
nextflow run isugifNF/parallel --input "ls *fastq | paste - -" --script "gsnapScript.sh" --threads 8 -profile nova
TranscriptomicsWorkshop/severin/Maize/02_gsnap/gsnapScript.sh" --threads 8 -profile nova
N E X T F L O W  ~  version 20.07.1
Launching `isugifNF/parallel` [scruffy_hoover] - revision: 101d292e2e [master]
process finished for inFILE
executor >  local (1), slurm (24)
[3c/8064fe] process > createInput      [100%] 1 of 1 ✔
[5a/10a85b] process > inputScript (17) [100%] 24 of 24 ✔
Completed at: 17-Feb-2021 14:03:30
Duration    : 2h 39m 27s
CPU hours   : 10.4
Succeeded   : 25
```
