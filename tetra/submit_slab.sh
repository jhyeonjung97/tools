#!/bin/bash
#SBATCH -N 2
#SBATCH -G 8
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -J jobname
#SBATCH -t 02:00:00
#SBATCH -A m2997
#SBATCH -e err.%j.log
#SBATCH -o out.%j.log

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load vasp/6.4.3-gpu

export VASP_SCRIPT=/global/homes/j/jiuy97/bin/run_vasp_gpu2.py
export VASP_PP_PATH=/global/cfs/cdirs/m2997/vasp-psp/pseudo54

python ~/bin/tools/tetra/opt_slab2_afm.py
python ~/bin/get_restart3.py