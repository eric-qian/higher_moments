#!/bin/bash
#SBATCH --job-name=mainSimul_T=500      # create a short name for your job
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=30        # cpu-cores per task
#SBATCH --mem-per-cpu=2G          # memory per cpu-core
#SBATCH --time=96:00:00           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all           # send email on job start, end and fault
#SBATCH --mail-user=eyqian@princeton.edu

module purge
module load matlab/R2019b

matlab -nodisplay -nosplash -r mainSimul