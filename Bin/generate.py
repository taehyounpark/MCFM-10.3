import os
import textwrap

import subprocess
from multiprocessing import Pool

sm_processes = ['gghZZ', 'gghZZ_x_ggZZ', 'ggZZ_box']
c6_processes = ['gghZZ', 'gghZZ_x_ggZZ']
c6_values = [-10.0, -5.0, -2.0, -1.0, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

def run_job(job):
    runstring, command = job
    # Define the log file path
    os.makedirs(runstring, exist_ok=True)
    logfilepath = os.path.join(runstring, f'{runstring}.out')
    """Execute a shell command."""
    # Open the file for writing and redirect stdout and stderr
    with open(logfilepath, 'w') as logfile:
        try:
            subprocess.run(command, shell=True, stdout=logfile, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed: {command}\nError: {e}")

def write_script(job):
    runstring, command = job
    script_contents = f"""#!/usr/bin/env bash
#SBATCH --job-name={runstring}
#SBATCH --output={runstring}/%j.out
#SBATCH --error={runstring}/%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G

{command}
"""
    script_path = f"{runstring}.sh"
    with open(script_path, 'w') as script_file:
        script_file.write(script_contents)
    return script_path

def submit_job(job):
    runstring, command = job
    os.makedirs(runstring, exist_ok=True)
    script_path = write_script(runstring, command)
    
    try:
        subprocess.run(f"sbatch {script_path}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit job {runstring}: {e}")

def define_jobs():
    """Generate a list of commands with different parameters."""
    jobs = []
    mcfm = './mcfm'

    for process in sm_processes:

        # sm
        runstring = f"{process}_sm"
        configfile = f"./input_{process}.ini"

        command = (
            f"{mcfm} {configfile} "
            f"-general%rundir={runstring} "
            f"-general%runstring={runstring} "
            f"-higgs_trilinear%c6=0.0"
        )
        jobs.append((runstring,command))

        # c6 values
    for process in c6_processes:
        for c6 in c6_values:
            runstring = f"{process}_c6_{str(c6).replace('.', 'd').replace('-', 'm')}"
            command = (
                f"{mcfm} {configfile} "
                f"-general%rundir={runstring} "
                f"-general%runstring={runstring} "
                f"-higgs_trilinear%c6={c6}"
            )
            jobs.append((runstring,command))
    
    return jobs

def main():
    jobs = define_jobs()
    # with Pool() as pool:
    #     pool.map(run_job, jobs)
    for job in jobs:
        submit_job(job)

if __name__ == '__main__':
    main()
