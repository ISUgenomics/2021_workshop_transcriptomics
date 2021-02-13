#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=gsnap
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=someone@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=project_name

set -e
set -u

start=`date +%s`

# === Load Modules here
# module load <name of module>

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Main Program

# Maize Reference Genome is B73... you may need to edit the gmap_build path
~/mybin/gmap_build -d b73 data_maize/ref/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna

end=`date +%s`

# === Log msgs and resource use                          
scontrol show job ${SLURM_JOB_ID}
echo "ran gsnap.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
