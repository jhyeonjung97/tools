#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1
#SBATCH -J IrFe-vib1
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

for dir in 1_Ir_top/1_V_V/ 1_Ir_top/6_O_OOH/ 2_Ir_hol/3_V_OH/ 2_Ir_hol/5_O_OH/ 4_IrFe_top2/6_O_OOH/ 5_IrFe_top3/6_O_OOH/ 6_IrFe_hol1/6_O_OOH/
do
    cp $dir/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/
done

