#!/bin/sh
#PBS -N jobname
#PBS -l select=1:ncpus=64:mpiprocs=16:ompthreads=4
#PBS -l walltime=06:00:00
#PBS -q normal
#PBS -A vasp
#PBS -V

cd $PBS_O_WORKDIR

#OpenMP settings:
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export VASP_SCRIPT=/home01/x2755a09/bin/tools/kisti/run_vasp16.py
export VASP_PP_PATH=/home01/x2755a09/POTCAR/

python /scratch/x2755a09/5_V_bulk/scripts/static_bulk2.py
python ~/bin/get_restart3

~/bin/lobster-5.1.0/lobster-5.1.0
python ~/bin/aloha/cohp.py > icohp.txt
python ~/bin/aloha/cobi.py > icobi.txt