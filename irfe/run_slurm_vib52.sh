#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=g2
#SBATCH -J IrFe-vib0
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

directories=(
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/5_Fe/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/5_Fe/2_V_O
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/5_Fe/3_V_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/5_Fe/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/5_Fe/7_OH_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/6_Fe/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/6_Fe/2_V_O
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/6_Fe/3_V_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/6_Fe/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/6_Fe/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/6_Fe/7_OH_OH
)

for dir in "${directories[@]}"; do
    cp $dir/vib/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/vib/
done