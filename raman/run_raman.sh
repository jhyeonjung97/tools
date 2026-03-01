#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1
#SBATCH -J raman-test
#SBATCH --time=05-00:00
#SBATCH -o out.%N.%j.log
#SBATCH -e err.%N.%j.log

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

export VASP_RAMAN_RUN="mpiexec.hydra -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.4.2/vasp.6.4.2.avx51
2.std.x"
export VASP_RAMAN_PARAMS='01_01_2_0.01'

python ~/bin/vasp_raman.py > vasp_raman.out