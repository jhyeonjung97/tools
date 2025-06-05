#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=g2
#SBATCH -J IrFe-vib1
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

directories=(
    /home/hyeonjung/scratch/4_IrFe3/1_H/0_Ir/3_vac_top
    /home/hyeonjung/scratch/4_IrFe3/1_H/0_Ir/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/2_Fe/3_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/2_Fe/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/3_Co/3_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/3_Co/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/4_Ni/3_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/4_Ni/4_vac_hol
)

for dir in "${directories[@]}"; do
    cp $dir/vib/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/vib/
done