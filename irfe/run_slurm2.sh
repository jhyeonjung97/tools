#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=g3
#SBATCH -J IrFe-R1
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

for dir in 1_Ir_top 2_Ir_hol 3_IrFe_top1 4_IrFe_top2 5_IrFe_top3
do
    for subdir in 1_V_V 2_V_O 3_V_OH 4_O_O 5_O_OH 6_O_OOH
    do
        cp $dir/$subdir/vib/* .
        mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
        cp * $dir/$subdir/vib/
    done
done