#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J test
#SBATCH -t 12:00:00
#SBATCH -A m2997
#SBATCH -e err.%j.log
#SBATCH -o out.%j.log

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load vasp-tpc/6.4.2-cpu

export VASP_SCRIPT=/global/homes/j/jiuy97/bin/run_vasp_cpu.py
export VASP_PP_PATH=/global/cfs/cdirs/m2997/vasp-psp/pseudo54

python ~/bin/tools/cep/mno2-100.py
python ~/bin/verve/bader.py
python ~/bin/get_restart3