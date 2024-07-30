import os
import textwrap

import subprocess
from multiprocessing import Pool

# sm_processes = ['gghZZ', 'gghZZ_x_ggZZ', 'ggZZ_box']
# c6_processes = ['gghZZ', 'gghZZ_x_ggZZ']
sm_processes = ['gghZZ']
c6_processes = ['gghZZ']
c6_values = [-10.0,-1.0, 1.0, 10.0]

def write_job(job):
    runstring, command = job
    os.makedirs(runstring, exist_ok=True)

    script_contents = f"""#!/usr/bin/env bash
#SBATCH --job-name={runstring}
#SBATCH --output={runstring}/%j.out
#SBATCH --error={runstring}/%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1       
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

source ./setup.sh
{command}
"""

    script_path = f"{runstring}/job.sh"
    with open(script_path, 'w') as script_file:
        script_file.write(script_contents)
    return script_path

def submit_job(job):
    runstring, command = job
    script_path = write_job(job)
    
    try:
        subprocess.run(f"sbatch {script_path}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit job {runstring}: {e}")

def define_job(process,c6):
    mcfm = './mcfm'
    configfile = f"./input_{process}.ini"
    runstring = f"{process}_sm" if c6==0.0 else f"{process}_c6_{str(c6).replace('.', 'd').replace('-', 'm')}"
    command = (
        f"{mcfm} {configfile} "
        f"-general%rundir={runstring} "
        f"-general%runstring={runstring} "
        f"-higgs_trilinear%c6={c6}"
    )
    return (runstring, command)

def main():

    jobs = []
    for process in sm_processes:
        jobs.append(define_job(process,0.0))
    for process in c6_processes:
        for c6 in c6_values:
            jobs.append(define_job(process,c6))

    for job in jobs:
        submit_job(job)

if __name__ == '__main__':
    main()
