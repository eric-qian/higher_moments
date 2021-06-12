#!/bin/bash
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=20       # cpu-cores per task
#SBATCH --mem-per-cpu=4G         # memory per cpu-core
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=eyqian@princeton.edu

module purge
module load matlab/R2019b

srun matlab -nodisplay -nosplash -r run_sim_var