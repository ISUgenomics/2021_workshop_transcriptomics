# 0 Data

* Log onto Atlas HPC, ([Atlas website: https://www.hpc.msstate.edu/computing/atlas/](https://www.hpc.msstate.edu/computing/atlas/))

```
ssh username@atlas-login.hpc.msstate.edu
```

* Check what nodes are available

```
sinfo
```

* Log onto the dtn (data transfer node) node, faster speed at transferring data onto the HPC

```
ssh atlas-dtn
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

