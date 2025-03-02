#!/bin/bash
#SBATCH -J jobname
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -G 4
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -A m2997
#SBATCH -e err.%j.log
#SBATCH -o out.%j.log

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load vasp/6.4.2-gpu

export VASP_SCRIPT=/global/homes/j/jiuy97/bin/run_vasp_gpu.py
export VASP_PP_PATH=/global/cfs/cdirs/m2997/vasp-psp/pseudo54

python ~/bin/tools/cation/opt_cluster.py
python ~/bin/verve/bader.py
python ~/bin/get_restart3