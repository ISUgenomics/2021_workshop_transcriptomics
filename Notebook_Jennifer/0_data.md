# 0 Data

## Transfer Maize & Bumblebee data to Atlas HPC

1. Log onto Atlas HPC, ([Atlas website: https://www.hpc.msstate.edu/computing/atlas/](https://www.hpc.msstate.edu/computing/atlas/))

  ```
  ssh username@atlas-login.hpc.msstate.edu
  ```

2. Check what nodes are available

  ```
  sinfo

  #> PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST 
  #> atlas*       up 14-00:00:0      1    mix Atlas-0001 
  #> atlas*       up 14-00:00:0    227   idle Atlas-[0002-0228] 
  #> bigmem       up 14-00:00:0      8   idle Atlas-[0229-0236] 
  #> gpu          up 14-00:00:0      4   idle Atlas-[0237-0240] 
  #> service      up 14-00:00:0      2   idle Atlas-dtn-[1-2] 
  ```

3. Log onto the dtn (data transfer node) node, faster speed at transferring data onto the HPC. Notice how the last partition `service` in the last column, has `Atlas-dtn-[1-2]`. This is our DTN node, the rest are compute nodes.

  ```
  ssh atlas-dtn     #<= will prompt you for password again
  ```

4. Either fetch Maize Data/Bee Data using `wget` or wrap it into a Slurm Script. The following fetches one file... we have many files to transfer. Ergo, instaed of waiting, we created a slurm script.

  ```
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_1.fastq.gz
  ```

  **Header for Slurm Script**

  ```
  #! /usr/bin/env bash
  #SBATCH --nodes=1
  #SBATCH -p service                  # <= partition (p) is dtn node    
  #SBATCH --ntasks-per-node=2
  #SBATCH --time=24:00:00
  #SBATCH --job-name=Maize
  #SBATCH --output=R-%x.%J.out
  #SBATCH --error=R-%x.%J.err
  #SBATCH --mail-user=username@email.com
  #SBATCH --mail-type=begin
  #SBATCH --mail-type=end
  #SBATCH --account=projectname       # <= required on Atlas

  # Place wget commands here.
  ```

## Maize Data

Since we had several Maize read files which could be fetched via ftp, we created a slurm script `atlas_maize.slurm`.

<details><summary>View `bin/atlas_maizedata.slurm`</summary>

```
#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH -p service            #<= notice how this is the dtn node
#SBATCH --ntasks-per-node=2
#SBATCH --time=24:00:00
#SBATCH --job-name=Maize
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=username@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=projectname       #<= atlas requires this

set -e
set -u

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Main Program
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_2.fastq.gz
# ... other wget commands, shortened for ease of reading

```

</details>

Submit a slurm script using `sbatch`.
 
```
sbatch atlas_maizedata.slurm
```

To verify it's running, look at the slurm queue 

```
squeue

#> JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
#> 130138_[1-13]     atlas maker-it joshua.u PD       0:00      1 (Dependency) 
#> 130140     atlas    busco joshua.u PD       0:00      1 (Dependency) 
#> 130142     atlas    blast joshua.u PD       0:00      1 (Dependency) 
#> 130143     atlas  iprscan joshua.u PD       0:00      1 (Dependency) 
#> 130137     atlas maker-it joshua.u PD       0:00      1 (DependencyNeverSatisfied) 
#> 130139     atlas maker-it joshua.u PD       0:00      1 (Dependency) 
#> 130141     atlas   filter joshua.u PD       0:00      1 (Dependency) 
#> 130120     atlas     bash brian.na  R      29:55      1 Atlas-0001 
```

Which lists all jobs on the queue. It's fun to see what other people are running. Your user name should be there. If `squeue` provides a list too long to read, I pipe to `grep` out my username.

```
# Only see my jobs on the queue
squeue | grep jennifer
```

For the Reference, we fetched the B73 Maize reference genome.

* NCBI Entry for Maize - [https://www.ncbi.nlm.nih.gov/assembly/GCF_902167145.1/](https://www.ncbi.nlm.nih.gov/assembly/GCF_902167145.1/)

  Fetch the fna (fasta nucleotide) and the gff (general feature format) files.

<!-- Need to unzip the fasta file before `gmap_build` can use it, or use the `--gunzip`. -->

  ```
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
  ```

## Bumblebee

1. Download from Google Drive... downloaded as several zipped files

  * [https://drive.google.com/drive/folders/17Avy8yNMwI0p-k0wv72x0uxaMd6UnnqF](https://drive.google.com/drive/folders/17Avy8yNMwI0p-k0wv72x0uxaMd6UnnqF)
  * Unzip the files (generates `0 Raw #` folders containing `fastq` files)

2. Move reads into one folder
  
  ```
  mkdir raw
  mv */*.fastq raw/.
  ```
  
3. copy to Atlas

  ```
  rsync -azP raw atlas-dtn:inbox/bee/.
  ```
  
4. For the Reference, we fetched the Bombus impatiens reference genome, found by `kemcelroy`.

  * [https://www.ncbi.nlm.nih.gov/genome/3415?genome\_assembly\_id=468201](https://www.ncbi.nlm.nih.gov/genome/3415?genome_assembly_id=468201)
  
  ```
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_genomic.fna.gz
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/GCF_000188095.3_BIMP_2.2_genomic.gff.gz
  ```
  
## Data Files

* Maize 
  * RNAseq reads = `/project/isu_gif_vrsc/2021-RNA_Workshop/data_maize/reads/*{1,2}fastq.gz`
  * reference = `/project/isu_gif_vrsc/2021-RNA_Workshop/data_maize/ref/`
* bumblebee 
  * RNAseq reads = `/project/isu_gif_vrsc/2021-RNA_Workshop/data_bee/reads/*.fastq`
  * reference = `/project/isu_gif_vrsc/2021-RNA_Workshop/data_bee/ref/`


## Quality Check (QC)

... in progress ...

```
module load fastqc
multiqc
```
