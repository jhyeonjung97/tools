#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1
#SBATCH -J IrFe-OH-brg
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

for dir in 0_Ir/4_Ir_brg 1_Mn/4_Ir_brg 2_Fe/4_Ir_brg 3_Co/4_Ir_brg 4_Ni/4_Ir_brg 1_Mn/6_M_brg 2_Fe/6_M_brg 3_Co/6_M_brg 4_Ni/6_M_brg
do
    cp $dir/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir
done