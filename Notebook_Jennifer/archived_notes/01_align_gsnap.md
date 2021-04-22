# 1 Align

Authors: Mou & Jennifer

Log onto Atlas HPC and check if the following programs are installed

* **GSNAP** - for indexing the reference and alignment
* **samtools** - for `*.sam` to `*.bam` file conversion (smaller intermediate files)
* **featureCounts** - for getting counts. 

## Check Atlas HPC Modules

The usual way to check if a program is installed is using:

```
module spider <programname>

#aka
module spider gsnap   # gmap? similar sounding name?
```

However this checks for a program one at a time. I prefer listing All modules that are available on an HPC and then using `less` or `grep` to search for the program or similarly named program.

Check what modules are available

```
(module avail) &> module_list.txt
```

Let's look at (`less`) the text file `module_list.txt`.

```
less module_list.txt
```

Which gives us a list of modules

```
------------------------------------------- /apps/modulefiles/core -------------------------------------------
   advisor/2020.2     gnuplot/5.2.7          mpich/3.3.2              python/3.7.5
   antlr/2.7.7        go/1.13.5              munge/0.5.13             qt/5.12.1
   boost/1.70         grace/5.1.25           namd/2.14                r/4.0.2
   bowtie2/2.4.1      gsl/2.6                ncl/6.6.2                rdesktop/2.2.0
   cairo/1.17.2       hdf5/1.10.6            nco/4.9.3                root/6.22.00
   canu/2.1           htslib/1.10.2          netcdf/4.7.4             rstudio/1.3.1073
   chapel/1.22.1      hwloc/2.2.0            nwchem/7.0.0             seqkit/0.14.0
   cmake/3.18.1       impi/2020.2            openblas/0.3.10          singularity/3.5.2
   contrib/0.1        intel/2020.2           openjdk/14.0.2           singularity/3.7.1 (D)
   cuda/11.0.3        intelpython3/2020.2    openmpi/4.0.4            slurm/20.02.4
   eigen/3.3.7        jasper/1.900.1         ovito/3.2.0              snopt/7.7
   exonerate/2.2.0    jellyfish/2.3.0        p4vasp/0.3.30            sparsehash/2.0.4
   fastqc/0.11.9      julia/1.5.1            parallel/20200722        sqlite/3.32.3
   ferret/7.5.0       kentutils/407          paraview/5.7.0           szip/2.1.1
   ffmpeg/4.3.1       local/0.1              pdsh/2.34                tassel/5.2.64
   fftw/3.3.8         make/4.3               perl/5.32.0              udunits/2.2.26
   gcc/10.2.0         matlab/2019b           pgi/2019-19.9            vep/100.4
   gd/2.3.0           mauve/2.4.0            pgi/2020-20.4     (D)    vmd/1.9.3
   gdal/3.1.2         meryl/1.0              pigz/2.4                 vtune/2020.2
   gdb/10.1           mesa/20.1.6            plink2/2.3               w3lib/2.0.6
   gempak/7.5.1       minimap2/2.17          prodigal/2.6.3           xcrysden/1.5.60
   geos/3.8.1         minpath/1.4            proj/7.1.0               xpdf/4.02
   git/2.28.0         mkl/2020.2             pycharm/2020.2           zlib/1.2.11

  Where:
   D:  Default Module

Use "module spider" to find all possible modules and extensions.
Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".
```

If we want to use any of the programs in this module list, we simply have to call

```
module load <module name>

# You could load a particular version module
module load singularity/3.5.2      #<= recommended to keep version number, for analysis reproducibility

# Or take the default version
module load singularity            #<= will load 3.7.1, since that has a "D" = default next to it in list
```

I'm not seeing GSNAP, samtools or featureCounts on the list...hmm, but there's `singularity` so we could install it via a singularity image. I'm surprised there isn't `miniconda` so we could install the python libraries. Mou and I worked on installing miniconda and GSNAP.

Side Note: To request something to be installed on Atlas HPC, email: `help-usda@hpc.msstate.edu`. Several programs are available via `singularity/3.5.2` but not the latest `singularity/3.7.1` yet...not sure when this will be updated.

<details><summary><b>Nova HPC</b> contained gsnap, samtools, and featureCounts</summary>

* [Nova\_module\_list.txt](bin/nova_module_list.txt)

Since Nova had many more modules, I used `grep` to pull out lines that contained "gsnap", "samtools" or "subread"

```
grep -e "gsnap" -e "samtools" -e "subread" nova_module_list.txt
```

Which gave me: 

```
   gmap-gsnap/2017-06-16-4oy56bt                        py-fastaindex/0.11rc7-py2-ilmrpiz                              r-rcpparmadillo/0.8.100.1.0-py2-m4g5gmx
   gmap-gsnap/2018-03-25-qa3kh3t                        py-fastaindex/0.11rc7-py2-sj3lhkk                       (D)    r-rcpparmadillo/0.8.100.1.0-py2-r3.5-6vrcyra
   gmap-gsnap/2018-07-04-gtu46xu                        py-faststructure/1.0-py2-k4vsldn                               r-rcpparmadillo/0.8.100.1.0-py2-r3.5-65z3l3h               (D)
   gmap-gsnap/2019-05-12-zjqshxf                 (D)    py-funcsigs/0.4-py2-prasuhx                                    r-rcppblaze/0.2.2-py2-r3.4-xd6bcfz
   help2man/1.47.4-phopsy7                              py-networkx/2.1-py2-fbsf2d3                                    r-rsamtools/1.28.0-py2-r3.4-cuda9-openmpi3-7qj6enh
   help2man/1.47.8-7sce2nu                              py-networkx/2.1-py3-mfsvnsu                                    r-rsamtools/1.32.2-py2-r3.5-mpich-hq3t6jr                  (D)
   libxml2/2.9.9-oqe2ao3                                r-a4reporting/1.24.0-py2-r3.4-wwnambt                          samtools/1.6-lyscjka
   libxml2/2.9.10-i7eqked                        (D)    r-acepack/1.4.1-py2-r3.4-7r7e46i                               samtools/1.7-kglvk7q
   libxmu/1.1.2-weujutd                                 r-acepack/1.4.1-py2-r3.5-3gq5fnp                               samtools/1.8-r54nmop
   libxmu/1.1.2-y6lkbh2                                 r-acepack/1.4.1-py2-r3.5-zmdesfl                        (D)    samtools/1.9-k6deoga
   libxmu/1.1.2-6zjibzx                          (D)    r-ade4/1.7-16-py3-pfdviww                                      samtools/1.10-py3-xuj7ylj                                  (D)
   meson/0.55.1-py3-wqljrz5                      (D)    r-blob/1.1.0-py2-r3.5-s7q5xdw                                  subread/1.6.0-ak6vxhs
   allinea/19.0.3    (D)    cplex/12.8-py2                   gmap-gsnap-legacy/2018.07.04        libs/fftw/3.3.4              perf-reports/6.0          starccm/13.04.010
```

Ergo, we only had to add the following to the top of slurm scripts:

```
# === Load Nova Modules
module load gmap-gsnap
module load samtools
module load subread

# === Can check the help documentation for each tool
gmap_build --help > gmap_build_help.txt
samtools view --help > samtools_help.txt
featureCounts --help > featureCounts_help.txt
```

</details>

---

## Install Programs

Before we get into installing your own programs in the HPC. We need to touch on file organization on an HPC.

----

### Side Tangent: Setting up folders in a HPC

This is my general template for organizing folders and files on an HPC. Since the home folder has a ~5GB memory limit, I tend to create folders in my project directory (larger memory) and then softlink (`ln -s`) them to home folder.  

```
/project/projectname/
           |_ software/              #<= soft link to home (shared programs across a research group, local installs here)
           |_ Jennifer/
                |_ inbox/            #<= softlink to home
                |_ outbox/           #<= softlink to home
                |_ dot_files/
                    |_ .singularity/ #<= softlink to home
                    |_ .conda/       #<= softlink to home
                    |_ .nextflow/    #<= softlink to home
                    |_ R/            #<= softlink to home
                    |_ miniconda3/   #<= will eventually add this... do not add this yet, could also put this in software
```

* `inbox` and `outbox`  are nice b.c. I can transfer files with `scp mylocalfile.tar.gz jenchang@atlas:inbox/.`  and don't have to think about 5gb limit
* dotfiles (`.singularity`, `.conda`) are usually invisible folders that get large as you install conda packages, or singularity images. These can eat up your home folder ~5GB memory limit if they're not softlinked

---

Initially we tried to install the programs natively, but eventually switched to `miniconda`

<details><summary>Notes from a local <b>GSNAP</b> install - WORKED</summary>

## Local install of GSNAP

Initially, I was indecisive on compiling GSNAP from source or installing via miniconda. Compiling from source is always going to run faster (more natively) then installing via miniconda/singularity/other containerization program. Installing from a miniconda/singularity/etc will tend to be easier. In the end, we were able to get it installed from source.

* GMAP and GSNAP Source code is available at first link - [http://research-pub.gene.com/gmap/](http://research-pub.gene.com/gmap/)

```
# ==== Place in shared software folder
cd /project/project_name/software   

# ==== Fetch GSNAP Source code
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-12-17.tar.gz

# ==== Unzip the tar.gz compressed folder
tar -xvzf gmap-gsnap-2020-12-17.tar.gz

# ==== Look around, notice there's a README, with an INSTALL file
cd gmap-2020-12-17
ls
#> acinclude.m4  AUTHORS	 config       configure     COPYING  LICENSE	  Makefile.in  NOTICE  src    util
#> aclocal.m4    ChangeLog  config.site  configure.ac  INSTALL  Makefile.am  NEWS	       README  tests  VERSION
```

Any C program with a Makefile can usually be installed with `./configure` -> `make` -> `make install`. 

If we look at the `README` instructions, it shows the `configure` -> `make` -> `make install` method of building.  Here, we're going to modify the `./configure` slightly to place the `gmap` executable in the `software` folder.

``` 
# ==== Compile
./configure --prefix=/project/project_name/software/gmap-2020-12-17 --with-gmapdb=/project/project_name/software/gmapdb
make
make install

# ==== Look around to see what has changed
ls
#> acinclude.m4  ChangeLog    config.status  INSTALL      Makefile.in  src
#> aclocal.m4    config	   configure	  LICENSE      NEWS	    tests
#> AUTHORS       config.log   configure.ac   Makefile     NOTICE	    util
#> bin	      config.site  COPYING	  Makefile.am  README	    VERSION
```

Difficult to see in the list, but now there's a `bin` folder. Look inside that folder and you will see the `gsnap` executables.

```
ls bin/

#> atoiindex	  gmap.avx2	gsnapl			    iit_store
#> cmetindex	  gmap_build	gsnapl.avx2		    indexdb_cat
#> cpuid		  gmap_cat	gsnapl.nosimd		    md_coords
#> dbsnp_iit	  gmapindex	gsnap.nosimd		    psl_genes
#> ensembl_genes	  gmapl		gtf_genes		    psl_introns
#> fa_coords	  gmapl.avx2	gtf_introns		    psl_splicesites
#> get-genome	  gmapl.nosimd	gtf_splicesites		    sam_sort
#> gff3_genes	  gmap.nosimd	gtf_transcript_splicesites  snpindex
#> gff3_introns	  gmap_process	gvf_iit			    trindex
#> gff3_splicesites  gsnap		iit_dump		    vcf_iit
#> gmap		  gsnap.avx2	iit_get
```
----

</details>

<details><summary>Notes from local install of <b>featureCounts</b> - WORKED</summary>

## Local install of featureCounts

Heh, this one was missing too.

* Install Instructions - [http://bioinf.wehi.edu.au/subread-package/](http://bioinf.wehi.edu.au/subread-package/)

Thankfully this is an already compiled binary. Ergo:

1. Download the `subread-2.0.1-Linux-x86_64.tar.gz` from sourceforge - [https://sourceforge.net/projects/subread/files/subread-2.0.1/](https://sourceforge.net/projects/subread/files/subread-2.0.1/)
2. scp file to Atlas

  ```
  # from local laptop, move file to Atlas
  scp subread-2.0.1-Linux-x86_64.tar.gz atlas:inbox/.
  ```
  
3. Move to software folder

  ```
  # From Atlas HPC
  mv ~/inbox/subread-2.0.1-Linux-x86_64.tar.gz /project/project_name/software/.
  cd /project/project_name/software
  tar -xzvf subread-2.0.1-Linux-x86_64.tar.gz
  cd subread-2.0.1-Linux-x86_64
  ls bin/*
  
  #> bin/exactSNP  bin/featureCounts  bin/subindel  bin/subjunc  bin/sublong  bin/subread-align  bin/subread-buildindex
  #> bin/utilities:
  #> detectionCall  flattenGTF  genRandomReads  propmapped  qualityScores  removeDup  repair  subread-fullscan  txUnique
  ```
  
4. Make sure program runs: `./bin/featureCounts`

</details>

<details><summary>Notes from a local install of **samtools** - c library linking errors, did not finish</summary>

## Local install of samtools

Geh... samtools is missing, which requires two difficult to install libraries... the benefits of installing natively has gone down. Decided to switch to miniconda. Samtools required a few specific C libraries, ergo I gave up natively install and switched to miniconda.

</details>


# Restart: install everything in miniconda

## Install miniconda

Based on [https://conda.io/projects/conda/en/latest/user-guide/install/linux.html](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html), install locally to Atlas HPC.

```
# Fetch the install script
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the fetched script
bash Miniconda3-latest-Linux-x86_64.sh
# acceept all defaults (press enter/yes through entire download process)

# By default it places miniconda3 in home folder (remember 5gb limit)
# Ergo we move this to the project folder and link to home
mv ~/miniconda3 /project/project_folder/software/.
ln -s /project/project_folder/software/miniconda3 ~/.

# Check if miniconda is installed and working
source ~/.bashrc
conda -version
```

Instead of moving the `miniconda3` folder, there are options of using the `--prefix` and `source activate path`.

* For more info - https://stackoverflow.com/questions/46929791/activating-conda-environment-with-its-full-path

## Install samtools, gmap, subread (featureCounts)

1. Create an environment file, save it as `gsnap_env.yml`.

  ```
  name: gsnap_env
  channels:
    - conda-forge
    - bioconda
    - defaults
  dependencies:
    - python=3.8
    - gmap 
    - samtools
    - subread
  ```

2. From the command line create the miniconda environment:

  ```
  conda env create -f gsnap_env.yml
  ```

  And enjoy the ease with which miniconda installs all dependencies for gsnap (`gmap`), samtools (`samtools`) and featureCounts (`featureCounts`). 
  
3. Check if environment was created:

  ```
  conda env list          #<= list all environments
  
  #> conda environments:
  #> base                  *  /home/jennifer.chang/miniconda3
  #> gsnap_env                /home/jennifer.chang/miniconda3/envs/gsnap_env
  ```
  
  Notice how a new folder is created in `miniconda` which contains all the python executables and python packages.
  
  ```
  ls /home/jennifer.chang/miniconda3/envs/gsnap_env/
  ```
  
4. Then the tools are available when you activate the environment.

  ```
  conda activate gsnap_env
  samtools --version        # check samtools version...
  gmap --version
  featureCounts -v
  ```

5. In order to use a miniconda environment in a slurm script, you'll need to activate yoru local miniconda, then the `gsnap_env`. Place the following after your `SBATCH` lines and before the `gsnap` / `samtools` / `featureCounts` calls.

  The top of the slurm script will require a few extra lines before `conda activate` will work.

  ```
  set +eu
  USER=jennifer.chang
  source /home/${USER}/miniconda3/etc/profile.d/conda.sh
  conda activate gsnap_env

  ## gsnap/samtools/featureCount commands here
  ```

-----

# Start GSNAP alignment

Hahaha, well we've finally made it (after a week of trying installation methods). 

Mou found a tutorial link on GSNAP:

* [https://hbctraining.github.io/Intro-to-rnaseq-hpc-gt/lessons/08\_rnaseq\_workflow.html](https://hbctraining.github.io/Intro-to-rnaseq-hpc-gt/lessons/08_rnaseq_workflow.html)

The general format of the pipeline is:

  ```
  (1) index genome -> (2) map reads to genome -> (3) get counts
  ```

The counts will be sent to Diffential Expression analysis programs.

## (1) Build reference genome index with `gmap_build`

1. Index the reference genome

<details><summary>More info on Maize Reference</summary>

Fetch Maize Reference (B73)

* NCBI Entry for Maize - [https://www.ncbi.nlm.nih.gov/assembly/GCF_902167145.1/](https://www.ncbi.nlm.nih.gov/assembly/GCF_902167145.1/)

Fetch the fna (fasta nucleotide) and the gff (general feature format) files.

Need to unzip the fasta file before `gmap_build` can use it. 

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
```

</details>

```
# ==== Activate miniconda
set +eu
USER=jennifer.chang
source /home/${USER}/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env

# ==== Connect the executable (either local or miniconda)
# GMAP_BUILD=/project/project_name/software/gmap-2020-12-17/bin/gmap_build
GMAP_BUILD=gmap_build

# ==== Define input/output variables
GENOME_NAME=b73
GENOME_FASTA=data_maize/ref/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
GMAPDB=/project/project_name/software/gmapdb

# ==== Main Run
${GMAP_BUILD} --gunzip -d ${GENOME_NAME} -D ${GMAPDB} ${GENOME_FASTA}
```

* [Maize_Runner.slurm](bin/Maize_Runner.slurm) - has been combined into one slurm script

Notice how we can split the `gmap_build` command in 3 sections. 

|flag | value | reason |
|:-:|:-:|:-:|
| -d | b73 | arbitrary name of the reference genome |
| -D | /project/project_name/software/gmapdb | directory where indexed genome will be saved |
|--gunzip| | input file is a compressed gz file |

How might you change the above command to run bumblebee data?

# (2) map RNAseq reads to reference with `gsnap`

```
#SBATCH --mem 300000

# ==== Activate miniconda
set +eu
USER=jennifer.chang
source /home/${USER}/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env

# ==== Link executables (local or miniconda)
# GSNAP=/project/project_name/software/gmap-2020-12-17/bin/gsnap
GSNAP=gsnap
SAMTOOLS=samtools

# ==== Define input/output variables
GENOME_NAME=b73
GMAPDB=/project/project_name/software/gmapdb
READNAME="Treatment_r1"
READONE=  # Path to left reads
READTWO=  # Path to right reads
OUTBAM=${READNAME}.Aligned.sortedByCoord.out.bam

# ==== Main Run
${GSNAP} \
 -d ${GENOME_NAME} -D ${GMAPDB} \
 -t 16 -M 2 -n 10 -N 1 \
 --quality-protocol=sanger -w 200000 --pairmax-rna=200000 \
 -E 1 -B 2 \
 -A sam \
 ${READONE} ${READTWO} | \
 ${SAMTOOLS} view -bS - > ${OUTBAM}
 
# ${SAMTOOLS} sort -m 3G - > ${OUTBAM} <= The sort statement keeps running out of memory and quitting, commenting out for now
```

**gsnap parameters**

|flag | value | reason |
|:-:|:-:|:-:|
| -d | b73 | arbitrary name of the reference genome |
| -D | /project/project_name/software/gmapdb | directory where indexed genome will be saved |
| -t | 16 | number of threads to use |
| -M | 2 | |
| -n | 10 | |
| -N | 1 | |
|--quality-protocol | sanger | use sanger (vs illumina) quality scores |
| -w | 200000 | window? |
| --pairmax-rna | 200000 | |
| -E | 1 | | 
| -B | 2 | |
| -A | sam | output a sam file? |

**samtools view parameters**

|flag | value | reason |
|:-:|:-:|:-:|
| -b | | output file as a bam file |
| -S | | auto detect if input is bam or sam |

The purpose of `samtools view -bS` is basically to convert a large sam file into a smaller bam file (binary)

* see explaination on SAM vs BAM here - [wikilink](https://en.wikipedia.org/wiki/SAMtools#:~:text=SAM%20files%20are%20human%2Dreadable,to%20work%20with%20than%20SAM.) - lol, whoever is updating the wiki is also using `[todo: ...]` statements :)

[todo: add shortened view of SAM file or example SAM file here ]

To determine amount of memory per node, use `sinfo`.

```
sinfo -O partition,allocnodes,memory
#         |          |           |_ size of memory per node in megabytes
#         |          |_ allowed allocating nodes
#         |_ name of partition (#SBATCH -p partitionNameHere)
```

For more information on `sinfo`, see [this link](http://manpages.ubuntu.com/manpages/cosmic/man1/sinfo.1.html)

### In progress, preparing input files

... In progress... preparing input files ... wasn't sure if `gsnap` automatically detects paired read data if I pass in a glob `*_1.fastq.gz` and `*_2.fastq.gz` or if it would do an all `_1` vs all `_2` comparison. Check this first before scaling up. I could hard code the command (make sure it's not being inefficient) but the developers of `gsnap` hopefully have fixed this potential issue in the design of the software (test and check).

Paired end reads must be fed into gsnap in order, basically print out the paired end reads.

```
ls data_maize/reads/* |\
  tr '\n' '\t'|\
  sed 's/_2.fastq.gz/_2.fastq.gz|/g'|\
  tr '|' '\n'|\
  awk '{print $1,$2}' > input.txt
```

Let's look at input.txt

<details><summary>input.txt</summary>

```
data_maize/reads/SRR1573504_1.fastq.gz data_maize/reads/SRR1573504_2.fastq.gz
data_maize/reads/SRR1573505_1.fastq.gz data_maize/reads/SRR1573505_2.fastq.gz
data_maize/reads/SRR1573506_1.fastq.gz data_maize/reads/SRR1573506_2.fastq.gz
data_maize/reads/SRR1573507_1.fastq.gz data_maize/reads/SRR1573507_2.fastq.gz
data_maize/reads/SRR1573508_1.fastq.gz data_maize/reads/SRR1573508_2.fastq.gz
data_maize/reads/SRR1573509_1.fastq.gz data_maize/reads/SRR1573509_2.fastq.gz
data_maize/reads/SRR1573510_1.fastq.gz data_maize/reads/SRR1573510_2.fastq.gz
data_maize/reads/SRR1573511_1.fastq.gz data_maize/reads/SRR1573511_2.fastq.gz
data_maize/reads/SRR1573512_1.fastq.gz data_maize/reads/SRR1573512_2.fastq.gz
data_maize/reads/SRR1573513_1.fastq.gz data_maize/reads/SRR1573513_2.fastq.gz
data_maize/reads/SRR1573514_1.fastq.gz data_maize/reads/SRR1573514_2.fastq.gz
data_maize/reads/SRR1573515_1.fastq.gz data_maize/reads/SRR1573515_2.fastq.gz
data_maize/reads/SRR1573516_1.fastq.gz data_maize/reads/SRR1573516_2.fastq.gz
data_maize/reads/SRR1573517_1.fastq.gz data_maize/reads/SRR1573517_2.fastq.gz
data_maize/reads/SRR1573518_1.fastq.gz data_maize/reads/SRR1573518_2.fastq.gz
data_maize/reads/SRR1573519_1.fastq.gz data_maize/reads/SRR1573519_2.fastq.gz
data_maize/reads/SRR1573520_1.fastq.gz data_maize/reads/SRR1573520_2.fastq.gz
data_maize/reads/SRR1573521_1.fastq.gz data_maize/reads/SRR1573521_2.fastq.gz
data_maize/reads/SRR1573522_1.fastq.gz data_maize/reads/SRR1573522_2.fastq.gz
data_maize/reads/SRR1573523_1.fastq.gz data_maize/reads/SRR1573523_2.fastq.gz
data_maize/reads/SRR1573524_1.fastq.gz data_maize/reads/SRR1573524_2.fastq.gz
data_maize/reads/SRR1573525_1.fastq.gz data_maize/reads/SRR1573525_2.fastq.gz
data_maize/reads/SRR1573526_1.fastq.gz data_maize/reads/SRR1573526_2.fastq.gz
data_maize/reads/SRR1573527_1.fastq.gz data_maize/reads/SRR1573527_2.fastq.gz
```

</details>

# (3) get counts with `featureCounts`

Posted by Rick in Slack Channel, will need to modify for gsnap output.

```
featureCounts -T 16 -p -t gene -g ID -a augustus.gff3 -o filtered_1703-TM102_sorted_counts_genes.txt filtered_1703-TM102_sorted.bam
```

**featureCounts parameters**

|flag | value | reason |
|:-:|:-:|:-:|
| -T | | threads?|
| -p | | |
| -t | gene| |
| -g | | |
| -a | augustus.gff3 | annotation?|
| -o | | output file name?|

## Run the above 3 steps with Bee data

* [Bee_Runner.slurm](bin/Bee_Runner.slurm)

Followed Mou's method of looping across files, ran on Atlas HPC. The next step is splitting this across several slurm jobs to run in parallel. The final output are named similar to `readname_genecounts.txt` which can either be combined in bash or in R as input to DE analysis programs.

```
#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Bees
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=jenchang@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=isu_gif_vrsc

set -e
set -u

start=`date +%s`

# === Load Modules here and link executables

# = Nova HPC
# module load gmap-gsnap
# module load samtools
# module load subread

# = Atlas HPC
set +eu
source /home/jennifer.chang/miniconda3/etc/profile.d/conda.sh
conda activate gsnap_env
GMAP_BUILD=gmap_build
GSNAP=gsnap
SAMTOOLS=samtools
FEATURECOUNTS=featureCounts

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Input / Output Variables
REF_NAME=Bombus
REF_FILE=data_bee/ref/GCF_000188095.3_BIMP_2.2_genomic.fna.gz
REF_GFF=data_bee/ref/GCF_000188095.3_BIMP_2.2_genomic.gff.gz
GMAPDB=gmapdb
# See forloop for the directory of reads

# === Main Program
# (1) Index Genome
${GMAP_BUILD} \
  --gunzip \
  -d ${REF_NAME} \
  -D ${GMAPDB} \
  ${REF_FILE}

for FILE in data_bee/reads/*.fastq
do
  READ_NAME=$(basename ${FILE} | sed 's:_L002_R1_001.fastq::g')
  DIR_NAME=$(dirname ${FILE})
  READ_R1=${DIR_NAME}/${READ_NAME}_L002_R1_001.fastq
  OUT_BAM=${READ_NAME}.aligned.out.bam
  OUT_COUNTS=${READ_NAME}_genecounts.txt

# (2) Map Reads:

  ${GSNAP} \
    -d ${REF_NAME} \
    -D ${GMAPDB} \
    -N 1 -t 16 -B 4 -m 5 \
    --input-buffer-size=1000000 \
    --output-buffer-size=1000000 \
    -A sam \
    ${READ_R1} | \
    ${SAMTOOLS} view --threads 16 -bS - > ${OUT_BAM}

# (3) Gene Counts
  ${FEATURECOUNTS} -T 16 -t gene -g ID \
    -a ${REF_GFF} \
    -o ${OUT_COUNTS} \
    ${OUT_BAM} 

done

end=`date +%s`

# === Log msgs and resource use                          
scontrol show job ${SLURM_JOB_ID}
echo "ran gsnap.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
```

# Counts

todo: describe final output here. Report basic stats, number of rows, etc. Do certain read pairs have more rows than others? Anything concerning about the data? Etc, etc.

```
# Program:featureCounts v2.0.1; Command:"featureCounts" "-T" "16" "-t" "gene" "-g" "ID" "-a" "data_bee/ref/GCF_000188095.3_BIMP_2.2_genomic.gff.gz" "-o" "1-A01-A1_S7_genecounts.txt" "1-A01-A1_S7.aligned.out.bam" 
Geneid  Chr     Start   End     Strand  Length  1-A01-A1_S7.aligned.out.bam
gene-LOC100740276       NT_176423.1     7       2256    +       2250    14
gene-LOC100740157       NT_176423.1     2829    5996    +       3168    100
gene-LOC100742884       NT_176427.1     27729   30739   +       3011    186
gene-LOC100740399       NT_176427.1     32165   37261   +       5097    25
gene-LOC100740519       NT_176427.1     38806   42290   -       3485    139
gene-LOC100743001       NT_176427.1     42433   53365   +       10933   112
gene-LOC100740639       NT_176427.1     54201   58114   +       3914    85
gene-LOC100743123       NT_176427.1     58465   60894   -       2430    149
...
```

Notice how the counts are in the final column.

Can combine in R.

* [combine.R](bin/combine.R)

```
#! /usr/bin/env Rscript
# Auth: Jennifer Chang
# Date: 2021/02/23
# Desc: Combine featureCounts output (1st and last column) files

# === Load Libraries
library(tidyverse)
library(magrittr)
library(readxl)

# === Get list of featureCount output files
dir_org="bee"         # counts are in a "bee" or "maize" subdirectory
featureCount_files <- list.files(path = dir_org, pattern = "*genecounts.txt$", full.names = TRUE)

# === Read in 1st file
data <- read_delim(featureCount_files[1], delim="\t", comment = "#" ) %>%
  select(Geneid, ends_with(".bam")) %>%              # Get 1st and last column (column was named after bam file)
  pivot_longer(cols=ends_with(".bam")) %>%           # Melt data (tidy data)
  mutate(
    name = gsub(".aligned.out.bam", "", name)        # No longer need the bam extension, easier to read
  )

# === Loop and append the rest
for (count_file in featureCount_files[-1]){
  print(count_file)
  temp <- read_delim(count_file, delim="\t", comment = "#") %>%
    select(Geneid, ends_with(".bam")) %>%
    pivot_longer(cols=ends_with(".bam")) %>%
    mutate(
      name = gsub(".aligned.out.bam", "", name)
    )
  data = rbind(data, temp)
}

# === Convert to excell like data (wider)
wide_data <- data %>%
  pivot_wider(id_cols=Geneid)

# === Save tab delimited file (smaller file size)
write_delim(wide_data, 
            paste(dir_org, "genecounts.txt", sep="_"), 
            delim="\t")

# === Save Excel file (can be easier to work with)
writexl::write_xlsx(wide_data, 
                    path=paste(dir_org, "genecounts.xlsx", sep="_"))
```

Final count files are in the following

* [bee_genecounts.txt](bee_genecounts.txt)
* [maize_genecounts.txt](maize_genecounts.txt)

The above two files can be fed into DESeq2 or EdgeR.
