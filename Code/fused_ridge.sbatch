#!/bin/sh

#SBATCH --output=slurmlogs/estimation_%x_%j.out
#SBATCH --error=slurmlogs/estimation_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5g
#SBATCH --mail-type=FAIL
#SBATCH --partition=broadwl


module load julia/1.6.2
julia fused_ridge_estimation.jl $SLURM_ARRAY_TASK_ID --thread 4
