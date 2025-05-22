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
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/1_layer_top
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/3_atom_top1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/5_atom_brg1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/6_atom_brg2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/6_IrFe/5_atom_brg1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/4_atom_top2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/5_atom_brg1
    /home/hyeonjung/scratch/4_IrFe3/2_OH/7_IrCo/6_atom_brg2
    /home/hyeonjung/scratch/4_IrFe3/2_OH/8_IrNi/5_atom_brg1
)

for dir in "${directories[@]}"; do
    cp $dir/vib/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * $dir/vib/
done