#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -q regular
#SBATCH -J jobname
#SBATCH -t 12:00:00
#SBATCH -A m2997
#SBATCH -e err.%j.log
#SBATCH -o out.%j.log

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load vasp-tpc/6.3.2-gpu

echo "import os" > run_vasp.py
echo "exitcode = os.system('srun -n 4 -c 32 vasp_std')" >> run_vasp.py

export VASP_SCRIPT=./run_vasp.py
export VASP_PP_PATH=/global/cfs/cdirs/m2997/vasp-psp/pseudo54

python /pscratch/sd/j/jiuy97/3_V_shape/scripts/static_bulk2.py
python ~/bin/verve/bader.py
python ~/bin/get_restart3

~/bin/lobster-5.0.0/lobster-5.0.0
python ~/bin/playground/aloha/cohp.py > icohp.txt
python ~/bin/playground/aloha/cobi.py > icobi.txt