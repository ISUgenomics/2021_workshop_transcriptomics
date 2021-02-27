#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=QC
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=jenchang@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e
set -u

start=`date +%s`

# === Load Modules here 
# source ~/mybin/config_module_names.sh        # <= switching between ceres and condo, may need to edit path

# = Nova HPC
module load parallel/20170322-36gxsog
module load fastqc/0.11.7-3flwcvl
module load py-multiqc/1.5-py2-lqqx3ht

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Input / Output
INDIR=Maize
OUTDIR=01_QC_${INDIR}

# === Main Program
[[ -d ${OUTDIR} ]] || mkdir ${OUTDIR}
parallel -j16 "fastqc -o ${OUTDIR}" ::: ${INDIR}/*.fastq
cd ${OUTDIR}
multiqc .

cd ..

# === Input / Output
INDIR=bumbleBee
OUTDIR=01_QC_${INDIR}

# === Main Program
[[ -d ${OUTDIR} ]] || mkdir ${OUTDIR}
parallel -j16 "fastqc -o ${OUTDIR}" ::: ${INDIR}/*.fastq
cd ${OUTDIR}
multiqc .

end=`date +%s`

# === Log msgs and resource use                          
scontrol show job ${SLURM_JOB_ID}
echo "ran QC.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
