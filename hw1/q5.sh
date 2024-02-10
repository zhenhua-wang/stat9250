#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH --job-name=q5_remote           # name for the job
#SBATCH -N 1                           # number of nodes
#SBATCH --tasks-per-node=1            # number of cores
#SBATCH --mem=4G                       # total memory
#SBATCH --time 0-00:10                 # time limit in the form days-hours:minutes
#SBATCH --mail-user=zwhkv@umsystem.edu # email address for notifications
#SBATCH --mail-type=FAIL
#SBATCH --partition=general
#--------------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled

echo "### Starting at: $(date) ###"
module purge
module load openmpi
module load miniconda3
source activate py

## Run the python script
for n_workers in 1 4 8 16 32 64
do
    echo
    echo "### Starting with $n_workers worker ###"
    Rscript ~/Workspace/stat9250/hw1/q5_remote.r $n_workers
    echo "### Ending with $n_workers worker ###"
    echo
done
conda deactivate
echo "### Ending at: $(date) ###"
