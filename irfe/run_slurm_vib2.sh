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
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/3_layer_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/4_layer_top2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/5_layer_top3
    /home/hyeonjung/scratch/4_IrFe3/2_OH/0_Ir/6_atom_top
    /home/hyeonjung/scratch/4_IrFe3/2_OH/2_Fe/5_layer_top3
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/3_layer_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/3_Co/5_layer_top3
    /home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/3_layer_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/4_Ni/4_layer_top2  
)

for dir in "${directories[@]}"; do
    cp $dir/vib/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/vib/
done