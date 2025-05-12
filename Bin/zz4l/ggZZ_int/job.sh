#!/usr/bin/env bash
#SBATCH --job-name=zz4l/ggZZ_int
#SBATCH --output=zz4l/ggZZ_int/%j.out
#SBATCH --error=zz4l/ggZZ_int/%j.err
#SBATCH --time=22:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

module purge
module load gcc/14

export OMP_STACKSIZE=16000
export OMP_NUM_THREADS=1

./mcfm ./input_ggZZ_int.ini -general%rundir=zz4l/ggZZ_int 
