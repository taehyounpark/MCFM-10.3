#!/usr/bin/env bash
#SBATCH --job-name=qqZZ2e2m
#SBATCH --output=qqZZ2e2m/%j.out
#SBATCH --error=qqZZ2e2m/%j.err
#SBATCH --time=22:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

module purge
module load gcc/14

export OMP_STACKSIZE=16000
export OMP_NUM_THREADS=1

./mcfm ./input_qqZZ2e2m.ini -general%rundir=qqZZ2e2m 
