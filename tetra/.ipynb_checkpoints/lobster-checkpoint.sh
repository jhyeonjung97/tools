#!/bin/sh
#PBS -N jobname
#PBS -l select=1:ncpus=64:mpiprocs=32:ompthreads=2
#PBS -l walltime=48:00:00
#PBS -q normal
#PBS -A vasp
#PBS -V

cd $PBS_O_WORKDIR

#OpenMP settings:
export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

~/bin/lobster-5.1.0/lobster-5.1.0