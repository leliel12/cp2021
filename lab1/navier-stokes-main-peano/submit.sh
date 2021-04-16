#!/bin/bash

#SBATCH --job-name=navier-stokes

#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpaz@unc.edu.ar

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --exclusive

srun perf stat -r 5 -dd -o perfO$1 ./headless_peano  $1> outO$1
