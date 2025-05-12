import os
import textwrap

import subprocess
from multiprocessing import Pool

basedir = 'zz4l'

processes = []
processes += ['ggZZ_sbi']
processes += ['ggZZ_sig']
processes += ['ggZZ_int']
processes += ['ggZZ_bkg']
processes += ['qqZZ']
# processes += ['qqWW']
# processes += ['ppZZ']

def write_job(job):
    rundir, command = job
    os.makedirs(rundir, exist_ok=True)

    script_contents = f"""#!/usr/bin/env bash
#SBATCH --job-name={rundir}
#SBATCH --output={rundir}/%j.out
#SBATCH --error={rundir}/%j.err
#SBATCH --time=22:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

module purge
module load gcc/14

export OMP_STACKSIZE=16000
export OMP_NUM_THREADS=1

{command}
"""

    script_path = f"{rundir}/job.sh"
    with open(script_path, 'w') as script_file:
        script_file.write(script_contents)
    return script_path

def submit_job(job):
    script_path = write_job(job)
    try:
        subprocess.run(f"sbatch {script_path}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit {script_path}: {e}")

def define_job(process):
    mcfm = './mcfm'
    configfile = f"./input_{process}.ini"
    rundir = f"{basedir}/{process}"
    command = (
        f"{mcfm} {configfile} "
        f"-general%rundir={rundir} "
    )
    return (rundir, command)

def main():

    jobs = []
    for process in processes:
        jobs.append(define_job(process))

    for job in jobs:
        submit_job(job)

if __name__ == '__main__':
    main()
