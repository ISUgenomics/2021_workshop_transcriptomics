# Hisat2 - alignment

HiSat2 alignment followed by Stringtie/Ballgown. This is definitely not done! Writing up notes before I start running commands.

## Install HiSat2

Visit the HISAT2 download page, fetch the linux version.

* [http://daehwankimlab.github.io/hisat2/download/](http://daehwankimlab.github.io/hisat2/download/)

Maybe secure copy (`scp`) it over to Atlas HPC...

```
scp ~/Downloads/hisat2-2.2.1-Linux_x86_64.zip atlas:inbox/.
```

Heh, maybe I should add a discussion of `~/.ssh/config` files. It shortens the HPC address, easier to navigate to.

## Maize

Hmm, does HISAT have a maize or bee genome? I should be able to build it manually right?

1. Index Genome

```
# === Input / Output Variables
REF_FILE=data_maize/ref/*.fna.gz
REF_NAME=b73

# === Main Program
hisat2-build ${REF_FILE} ${REF_NAME}
```

Will generate genome index files with a `*.ht21` file extension.

2. Run HiSAT to get counts?

Hmm... not a fan of this input format... maybe there's a param I'm missing.

Still reading through the documentation...

* [http://daehwankimlab.github.io/hisat2/manual/](http://daehwankimlab.github.io/hisat2/manual/)

```
# === Input / Output Variables
REF_FILE=data_maize/ref/*.fna.gz
REF_NAME=b73

# === Main Program
# create comma separated list of left and right reads...
ls data_maize/reads/*_1.fastq.gz |tr '\n' ',' > reads_1.txt
ls data_maize/reads/*_2.fastq.gz |tr '\n' ',' > reads_2.txt

hisat2 -x ${REF_NAME} \
  -1 $(cat reads_1.txt) \
  -2 $(cat reads_2.txt) \
  -p 16 \
  -S ${READ_NAME}.sam
```

Can I redirect this to `samtools view` to get a smaller bam file? Wait is this combining all together? Maybe I still need the file loop.
