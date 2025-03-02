#!/bin/bash
#SBATCH -J pt
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 120
#SBATCH -p roma
#SBATCH --mem=480G
#SBATCH -e err.%j.log
#SBATCH -o out.%j.log
#SBATCH -A suncat:normal

module load mpi

export OMP_NUM_THREADS=1
export VASP_SCRIPT=/sdf/home/j/jiuy97/bin/run_vasp_cpu.py
export VASP_PP_PATH=/fs/ddn/sdf/group/suncat/sw/psp/vasp/potpaw54

python ./opt_slab.py
python ~/bin/verve/bader.py
python ~/bin/get_restart3