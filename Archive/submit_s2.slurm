#!/bin/bash
#SBATCH -J s2
#SBATCH -c 2
#SBATCH -t 30-0

## Sepcify the working directory for this job
cd $SLURM_SUBMIT_DIR

source /app/lmod/lmod/init/profile
module load fhR/4.4.0-foss-2023b

Rscript comboAMP_exponentialModel_iterated_s2.R

exit 0
