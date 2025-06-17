#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=g2
#SBATCH -J IrFe-vib
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

directories=(
    /home/hyeonjung/scratch/4_IrFe3/1_H/1_Mn/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/3_Fe/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/3_Fe/7_OH_OH
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/4_Fe/2_V_O
)

for dir in "${directories[@]}"; do
    cp $dir/vib/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/vib/
done