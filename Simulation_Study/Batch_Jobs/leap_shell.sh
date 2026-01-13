#!/bin/bash

#SBATCH --job-name=leap
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH	--array=36-40,56-60,151-153,154-160
#SBATCH --output=/work/users/c/l/clairez1/Paper1sims3/LEAP_model/Log/slurmLogFiles%a.out
#SBATCH --error=/work/users/c/l/clairez1/Paper1sims3/LEAP_model/Error/%a.err
#SBATCH --constraint=rhel8

## add R module
module add gcc/11.2.0
module add r/4.3.1

R CMD BATCH --no-restore /nas/longleaf/home/clairez1/Paper1_sims_rev2/simulation_programs/leap_sim.R /work/users/c/l/clairez1/Paper1sims3/LEAP_model/sim_$SLURM_ARRAY_TASK_ID.Rout


