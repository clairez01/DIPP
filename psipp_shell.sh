#!/bin/bash

#SBATCH --job-name=psipp_sims
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3g
#SBATCH	--array=1-600
#SBATCH --output=/work/users/c/l/clairez1/Paper1sims2/PSIPP/Log/slurmLogFiles%a.out
#SBATCH --error=/work/users/c/l/clairez1/Paper1sims2/PSIPP/Error/%a.err
#SBATCH --constraint=rhel8

## add R module
module add gcc/11.2.0
module add r/4.3.1

R CMD BATCH --no-restore /nas/longleaf/home/clairez1/R/Paper1sims2_files/nimble_psipp_simcode.R /work/users/c/l/clairez1/Paper1sims2/PSIPP/sim_$SLURM_ARRAY_TASK_ID.Rout


