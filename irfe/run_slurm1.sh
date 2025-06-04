#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1
<<<<<<< HEAD
#SBATCH -J IrFe-H1
=======
#SBATCH -J IrFe-OH
>>>>>>> 733760c112ad7de7bbd9703d39c54f2d49105166
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

<<<<<<< HEAD
for dir in 3_vac_hol 4_vac_hol
do
    cp ../$dir/* .
    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
    cp * ../$dir/
=======
for dir in 1_V_V 2_V_O 3_V_OH 4_O_O 5_O_OH 6_O_OOH 7_OH_OH
do
    if [[ -f $dir/POSCAR ]]; then
        cp $dir/* .
        mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
        cp * $dir/
    fi
>>>>>>> 733760c112ad7de7bbd9703d39c54f2d49105166
done
