#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=g3
#SBATCH -J IrFe-O
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

for dir in 5_IrMn 6_IrFe
do
    for subdir in 2_layer_hol 3_atom_top1 5_atom_hol1
    do
        cp $dir/$subdir/* .
        mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
        cp * $dir/$subdir/
    done
done

for dir in 7_IrCo 8_IrNi
do
    for subdir in 3_atom_top1 5_atom_hol1
    do
        cp $dir/$subdir/* .
        mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
        cp * $dir/$subdir/
    done
done