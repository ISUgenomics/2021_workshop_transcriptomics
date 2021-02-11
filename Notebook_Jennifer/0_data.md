# 0 Data

* Log onto Atlas HPC, ([Atlas website: https://www.hpc.msstate.edu/computing/atlas/](https://www.hpc.msstate.edu/computing/atlas/))

```
ssh username@atlas-login.hpc.msstate.edu
```

* Check what nodes are available

```
sinfo

#> PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST 
#> atlas*       up 14-00:00:0      1    mix Atlas-0001 
#> atlas*       up 14-00:00:0    227   idle Atlas-[0002-0228] 
#> bigmem       up 14-00:00:0      8   idle Atlas-[0229-0236] 
#> gpu          up 14-00:00:0      4   idle Atlas-[0237-0240] 
#> service      up 14-00:00:0      2   idle Atlas-dtn-[1-2] 
```

* Log onto the dtn (data transfer node) node, faster speed at transferring data onto the HPC. Notice how the last partition `service` in the last column, has `Atlas-dtn-[1-2]`. This is our DTN node, the rest are compute nodes.

```
ssh atlas-dtn     #<= will prompt you for password again
```

* Either fetch Maize Data/Bee Data using `wget` or wrap it into a Slurm Script:

```
sbatch atlas_maizedata.slurm
```

<details><summary>View **bin/atlas_maizedata.slurm**</summary>

```
#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH -p service         #<= notice how this is the dtn node
#SBATCH --ntasks-per-node=2
#SBATCH --time=24:00:00
#SBATCH --job-name=Maize
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
# --mail-user=username@email.com
# --mail-type=begin
# --mail-type=end
# --account=projectname   #<= atlas requires this

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

Which lists all jobs no the queue. It's fun to see what other people are running. Your user name should be there. If `squeue` provides a list too long to read, I grep out my username.

```
# Only see my jobs on the queue
squeue | grep jennifer
```
