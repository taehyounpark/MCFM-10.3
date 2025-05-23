import os
import textwrap

import subprocess
from multiprocessing import Pool

# mode = 'zz4l'
mode = 'zz2l2v'

processes = []
processes += ['ggZZ_sbi']
processes += ['ggZZ_sig']
processes += ['ggZZ_int']
processes += ['ggZZ_bkg']
processes += ['qqZZ']
processes += ['qqWW']
# processes += ['ppZZ']

def write_job(mode, proc):
    mcfm = './mcfm'
    
    rundir = os.path.join(mode,proc)
    os.makedirs(rundir, exist_ok=True)

    cfg = os.path.join(mode,f"input_{proc}.ini")
    cmd = f"{mcfm} {cfg} -general%runstring={mode} -general%rundir={rundir}"

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

    script_path = write_job(mode, proc)

    try:
        subprocess.run(f"sbatch {script_path}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit {script_path}: {e}")

def main():

    for process in processes:
        submit_job(mode, process)

if __name__ == '__main__':
    main()
