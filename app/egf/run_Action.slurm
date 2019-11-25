#!/bin/bash
#SBATCH -J Action
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -t 72:00:00
#SBATCH --qos=normal
#SBATCH -o slurm-%A.out


module add matlab
matlab -nojvm -nodesktop -nodisplay -nosplash -r 'Action'
