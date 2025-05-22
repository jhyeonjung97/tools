#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=g2
#SBATCH -J IrFe
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

directories=(
    /home/hyeonjung/scratch/4_IrFe3/4_R1/1_Ir_top/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/4_R1/2_Ir_hol/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/4_R1/2_Ir_hol/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/4_R1/3_IrFe_top1/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/4_R1/4_IrFe_top2/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/4_R1/4_IrFe_top2/6_O_OOH
)

for dir in "${directories[@]}"; do
    cp $dir/vib/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/vib/
done