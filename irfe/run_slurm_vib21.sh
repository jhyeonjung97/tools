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
)

for dir in "${directories[@]}"; do
    cp $dir/vib/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/vib/
done