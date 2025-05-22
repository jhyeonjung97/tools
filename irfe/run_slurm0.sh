#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1
#SBATCH -J IrFe
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

directories=(
    /home/hyeonjung/scratch/4_IrFe3/1_H/1_Mn/3_Ir_top
    /home/hyeonjung/scratch/4_IrFe3/1_H/3_Co/2_layer_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/4_Ni/2_layer_hol
    /home/hyeonjung/scratch/4_IrFe3/1_H/7_IrCo/1_layer_top
    /home/hyeonjung/scratch/4_IrFe3/2_OH/1_Mn/2_layer_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/1_Mn/4_Ir_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/2_layer_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/4_Ir_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/6_M_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/2_layer_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/3_Ir_top
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/4_Ir_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/6_M_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/2_layer_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/4_Ir_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/6_M_brg
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/1_layer_top
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/3_atom_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/5_atom_brg1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/6_atom_brg2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/6_IrFe/5_atom_brg1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/4_atom_top2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/5_atom_brg1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/6_atom_brg2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/8_IrNi/5_atom_brg1
    /home/hyeonjung/scratch/4_IrFe3/3_O/5_IrMn/3_atom_top1
    /home/hyeonjung/scratch/4_IrFe3/3_O/5_IrMn/5_atom_hol1
    /home/hyeonjung/scratch/4_IrFe3/3_O/6_IrFe/3_atom_top1
    /home/hyeonjung/scratch/4_IrFe3/3_O/7_IrCo/3_atom_top1
    /home/hyeonjung/scratch/4_IrFe3/3_O/8_IrNi/3_atom_top1
    /home/hyeonjung/scratch/4_IrFe3/3_O/8_IrNi/5_atom_hol1
    /home/hyeonjung/scratch/4_IrFe3/4_R1/1_Ir_top/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/4_R1/2_Ir_hol/5_O_OH
    /home/hyeonjung/scratch/4_IrFe3/4_R1/2_Ir_hol/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/4_R1/3_IrFe_top1/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/4_R1/4_IrFe_top2/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/4_R1/4_IrFe_top2/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/5_R2/1_Ir_top/1_V_V
    /home/hyeonjung/scratch/4_IrFe3/5_R2/2_Ir_hol/6_O_OOH
    /home/hyeonjung/scratch/4_IrFe3/5_R2/3_IrFe_top1/1_V_V
)

for dir in "${directories[@]}"; do
    cp $dir/CONTCAR ./POSCAR
    cp $dir/INCAR .
    cp $dir/KPOINTS .
    cp $dir/POTCAR .
    sed -i '/NSW/c\NSW = 200' INCAR
    touch continue.log
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/
done