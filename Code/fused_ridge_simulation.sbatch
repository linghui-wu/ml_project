#!/bin/sh

#SBATCH --output=slurmlogs/estimation_%x_%j.out
#SBATCH --error=slurmlogs/estimation_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5g
#SBATCH --mail-type=FAIL
#SBATCH --partition=broadwl
#SBATCH --array=1-10:1

module load julia/1.6.2
julia fused_ridge_simulation.jl const $SLURM_ARRAY_TASK_ID --thread 4
julia fused_ridge_simulation.jl independent $SLURM_ARRAY_TASK_ID --thread 4
