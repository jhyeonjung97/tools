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
    /home/hyeonjung/scratch/4_IrFe3/1_H/0_Ir/3_vac_top
    /home/hyeonjung/scratch/4_IrFe3/1_H/0_Ir/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/2_Fe/3_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/2_Fe/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/3_Co/3_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/3_Co/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/4_Ni/3_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/4_Ni/4_vac_hol
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/3_layer_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/4_layer_top2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/5_layer_top3
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/6_atom_top
    /home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/5_layer_top3
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/3_layer_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/5_layer_top3
    /home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/3_layer_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/4_layer_top2
    /home/hyeonjung/scratch/4_IrFe3/3_O/6_IrFe/2_layer_hol
    /home/hyeonjung/scratch/4_IrFe3/3_O/8_IrNi/2_layer_hol
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/2_Ir/7_OH_OH
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/3_Fe/2_V_O
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/3_Fe/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/4_Fe/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/4_Fe/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/5_Fe/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/5_Fe/2_V_O
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/5_Fe/7_OH_OH
    /home/hyeonjung/scratch/4_IrFe3/4_R_top/6_Fe/7_OH_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/3_Fe/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/3_Fe/2_V_O
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/3_Fe/3_V_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/3_Fe/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/3_Fe/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/3_Fe/7_OH_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/4_Fe/2_V_O
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/4_Fe/3_V_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/4_Fe/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/4_Fe/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/5_R_hol/4_Fe/7_OH_OH
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