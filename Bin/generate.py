import os
import textwrap

import subprocess
from multiprocessing import Pool

mode = 'zz4l'

processes = []
processes += ['ggZZ_sbi']
processes += ['ggZZ_sig']
processes += ['ggZZ_int']
processes += ['ggZZ_bkg']
processes += ['qqZZ']
# processes += ['qqWW']
# processes += ['ppZZ']

def form_command(mode, proc):
    mcfm = './mcfm'
    cfg = os.path.join(mode,f"input_{proc}.ini")
    command = f"{mcfm} {cfg} "
    return command

def write_script(mode, proc):
    cmd = form_command(mode, proc)

    rundir = os.path.join(mode,proc)
    os.makedirs(dir, exist_ok=True)

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

{cmd}
"""

    script_path = f"{rundir}/job.sh"
    with open(script_path, 'w') as script_file:
        script_file.write(script_contents)
    return script_path

def submit_job(mode, proc):

    script_path = write_script(mode, proc)

    try:
        subprocess.run(f"sbatch {script_path}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit {script_path}: {e}")

def main():

    for process in processes:
        submit_job(mode, process)

if __name__ == '__main__':
    main()
