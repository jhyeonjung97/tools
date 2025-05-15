#!/bin/bash

directories=(
    "/home/hyeonjung/scratch/4_IrFe3/1_H/1_Mn/3_Ir_top"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/3_Ir_top"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/1_Mn/3_Ir_top"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/1_Mn/4_Ir_brg"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/1_Mn/6_M_brg"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/2_layer_brg"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/3_Ir_top"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/4_Ir_brg"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/3_Ir_top"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/4_Ir_brg"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/3_Ir_top"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/4_atom_top2"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/6_atom_brg2"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/6_IrFe/5_atom_brg1"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/4_atom_top2"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/5_atom_brg1"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/6_atom_brg2"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/8_IrNi/5_atom_brg1"
    "/home/hyeonjung/scratch/4_IrFe3/2_OH/8_IrNi/6_atom_brg2"
    "/home/hyeonjung/scratch/4_IrFe3/3_O/1_Mn/3_Ir_top"
    "/home/hyeonjung/scratch/4_IrFe3/3_O/5_IrMn/3_atom_top1"
    "/home/hyeonjung/scratch/4_IrFe3/4_R1/1_Ir_top/1_V_V"
    "/home/hyeonjung/scratch/4_IrFe3/4_R1/2_Ir_hol/5_O_OH"
    "/home/hyeonjung/scratch/4_IrFe3/4_R1/2_Ir_hol/6_O_OOH"
    "/home/hyeonjung/scratch/4_IrFe3/4_R1/3_IrFe_top1/1_V_V"
    "/home/hyeonjung/scratch/4_IrFe3/4_R1/3_IrFe_top1/6_O_OOH"
    "/home/hyeonjung/scratch/4_IrFe3/4_R1/4_IrFe_top2/1_V_V"
    "/home/hyeonjung/scratch/4_IrFe3/4_R1/4_IrFe_top2/6_O_OOH"
    "/home/hyeonjung/scratch/4_IrFe3/5_R2/1_Ir_top/1_V_V"
    "/home/hyeonjung/scratch/4_IrFe3/5_R2/2_Ir_hol/6_O_OOH"
    "/home/hyeonjung/scratch/4_IrFe3/5_R2/3_IrFe_top1/1_V_V"
    "/home/hyeonjung/scratch/4_IrFe3/5_R2/3_IrFe_top1/6_O_OOH"
    "/home/hyeonjung/scratch/4_IrFe3/5_R2/5_IrFe_top3/6_O_OOH"
)

for dir in "${directories[@]}"; do
    cd "$dir" || continue
    if [[ -f OUTCAR ]] && [[ ! -f DONE ]]; then
        # touch continue.log
        # cp CONTCAR POSCAR
        pwd
    fi
done