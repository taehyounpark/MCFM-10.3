import os
import argparse
import subprocess

# Define your processes here
# mode = 'zz4l'

processes = {
    'zz4l' : [
        # 'ggZZ_sbi',
        # 'ggZZ_sig',
        'ggZZ_int',
        # 'ggZZ_bkg',
        # 'qqZZ'
        ],
    'zz2l2v' : [
        # 'ggZZ_sbi',
        # 'ggZZ_sig',
        'ggZZ_int',
        # 'ggZZ_bkg',
        # 'qqZZ'
        ],
}

def write_job(mode, proc, nthreads, time_hours, memory_gb):
    mcfm = './mcfm'
    rundir = os.path.join(mode, proc)
    os.makedirs(rundir, exist_ok=True)

    cfg = os.path.join(mode, f"input_{proc}.ini")
    cmd = f"{mcfm} {cfg} -general%runstring={mode} -general%rundir={rundir}"

    script_contents = f"""#!/usr/bin/env bash
#SBATCH --job-name={rundir}
#SBATCH --output={rundir}/%j.out
#SBATCH --error={rundir}/%j.err
#SBATCH --time={time_hours}:00:00
#SBATCH --mem={memory_gb}G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={nthreads}
#SBATCH --partition=general

module purge
module load gcc/14

export OMP_STACKSIZE=16000
export OMP_NUM_THREADS={nthreads}

{cmd}
"""

    script_path = f"{rundir}/job.sh"
    with open(script_path, 'w') as script_file:
        script_file.write(script_contents)
    return script_path

def submit_job(mode, proc, nthreads, time_hours, memory_gb):
    script_path = write_job(mode, proc, nthreads, time_hours, memory_gb)
    try:
        subprocess.run(f"sbatch {script_path}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit {script_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Submit MCFM jobs to SLURM.")
    parser.add_argument('mode', help='Run mode (e.g., zz2l2v, zz4l)')
    parser.add_argument('--nthreads', '-j', type=int, default=1, help='Number of OpenMP threads')
    parser.add_argument('--time', '-t', type=int, default=6, help='Walltime in hours')
    parser.add_argument('--memory', '-m', type=int, default=8, help='Memory in GB')
    args = parser.parse_args()

    for proc in processes[args.mode]:
        submit_job(args.mode, proc, args.nthreads, args.time, args.memory)

if __name__ == '__main__':
    main()
