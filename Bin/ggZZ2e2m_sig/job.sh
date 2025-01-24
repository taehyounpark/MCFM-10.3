#!/usr/bin/env bash
#SBATCH --job-name=ggZZ2e2m_sig
#SBATCH --output=ggZZ2e2m_sig/%j.out
#SBATCH --error=ggZZ2e2m_sig/%j.err
#SBATCH --time=22:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

module purge
module load gcc/14

export OMP_STACKSIZE=16000
export OMP_NUM_THREADS=1

./mcfm ./input_ggZZ2e2m_sig.ini -general%rundir=ggZZ2e2m_sig 
