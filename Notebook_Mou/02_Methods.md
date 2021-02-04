# Methods

Written summary of methods performed in this repo. This is the methods write up for the paper.

## HPC
Atlas dtn node
1. `ssh username@Atlas-dtn.hpc.msstate.edu`
https://www.hpc.msstate.edu/computing/atlas/
2. Project path:
* fsepru/username/rnaseq/maize_sequences
  * Directions on how to wget maize sequences: https://github.com/ISUgenomics/2021_workshop_transcriptomics/blob/main/00_Files.md
```
#!/bin/bash
#SBATCH --job-name=ncbi                              # name of the job submitted
#SBATCH -p service                                   # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 2                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 12:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#SBATCH --account fsepru
#Enter commands here:

#End of file
```
* To do: slurm script for downloading files

  * **change http to ftp**
* fsepru/username/rnaseq/bee_sequences


## Raw data
Maize: https://www.ebi.ac.uk/ena/browser/view/PRJNA260793
Bee data


## Sequence Analysis
1. Download Maize data
2. Download Bee data
