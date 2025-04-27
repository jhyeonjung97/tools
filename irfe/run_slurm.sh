#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1
#SBATCH -J IrFe-vib1
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

for dir in 5_IrMn 6_IrFe 7_IrCo 8_IrNi
do
    for subdir in 1_layer_top 2_layer_hol 3_atom_top1 4_atom_top2 5_atom_hol1 6_atom_hol2
    do
        if [[ -d ../$dir/$subdir/vib ]] && [[ ! -f ../$dir/$subdir/vib/OUTCAR ]]; then
            cp ../$dir/$subdir/vib/* .
            mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
            cp * ../$dir/$subdir/vib/
        fi
    done
done
