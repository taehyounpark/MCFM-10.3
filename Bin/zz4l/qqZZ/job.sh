#!/usr/bin/env bash
#SBATCH --job-name=zz4l/qqZZ
#SBATCH --output=zz4l/qqZZ/%j.out
#SBATCH --error=zz4l/qqZZ/%j.err
#SBATCH --time=22:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

module purge
module load gcc/14

export OMP_STACKSIZE=16000
export OMP_NUM_THREADS=1

./mcfm ./input_qqZZ.ini -general%rundir=zz4l/qqZZ 
