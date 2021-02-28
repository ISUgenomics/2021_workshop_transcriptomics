# run gsnap on Bumble Bee

* Feb 18, 2021
* /work/gif/TranscriptomicsWorkshop/severin/Bee/02_gsnap

## Softlink raw data and genome and gff

```
ls /work/gif/TranscriptomicsWorkshop/bumbleBee/*gz  | xargs -I xx ln -s xx
ln -s /work/gif/TranscriptomicsWorkshop/bumbleBee/01_Genome/GCF_000188095.3_BIMP_2.2_genomic.fna
cp /work/gif/TranscriptomicsWorkshop/bumbleBee/01_Genome/GCF_000188095.3_BIMP_2.2_genomic.gff.gz .
gunzip GCF_000188095.3_BIMP_2.2_genomic.gff.gz

```


## install latest gsnap/gmap
Modules are  out of date so installed it from scratch.  

```
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-12-17.tar.gz
tar -zxvf gmap-gsnap-2020-12-17.tar.gz
./configure --prefix=`pwd`
make
make install
# gmap_build kept looking for a bin directory so I created one and softlinked the executable files in util and src
mkdir bin; cd bin
ls -lha ../src | grep "r-xr-x" | xargs -I xx ln -s ./src/xx
ls -lha ../util | grep "r-xr-x" | xargs -I xx ln -s ./util/xx

```

### Build the Maize database for gsnap/gmap

* /work/gif/TranscriptomicsWorkshop/severin/Bee/02_gsnap

```
#module load gmap-gsnap  # version 2019-05-12-zjqshxf
salloc -N 1 -n 36 -p huge -t 4:00:00
export PATH=/work/gif/TranscriptomicsWorkshop/severin/gmap-2020-12-17/src/:$PATH
gmap_build -D BB -d BB GCF_000188095.3_BIMP_2.2_genomic.fna
```

### Run GSNAP using Nextflow

####  gsnapScript.sh
```
#!/bin/bash

module load gmap-gsnap/2018-07-04-gtu46xu
export GMAPDB=/work/gif/TranscriptomicsWorkshop/severin/Bee/02_gsnap/BB/BB/
DB_NAME="BB"
FILE1="$1"
FILE2="$2"
OUTFILE=$(basename ${FILE1} | sed 's/_1.fastq$//g')
# Note: "-N" option for detecting novel splice sites, remove if not needed (0=OFF; 1=ON)
gsnap -d ${DB_NAME} -N 1 -t 8 -B 4 -m 5 --input-buffer-size=1000000 --output-buffer-size=1000000 -A sam --split-output=${DB_NAME}_${OUTFILE} ${FILE1} ${FILE2}
```

#### run the script using nextflow parallel workflow

Parallel workflow requires the input files command which will be executed and placed in a queue that will be fed to the script parameter which in our case is the script above.

```
nextflow run isugifNF/parallel --input "ls /work/gif/TranscriptomicsWorkshop/severin/Bee/02_gsnap/*fastq.gz" --script "/work/gif/TranscriptomicsWorkshop/severin/Bee/02_gsnap/gsnapScript.sh" --threads 8 -profile nova

```
