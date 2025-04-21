#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1
#SBATCH -J IrFe1
#SBATCH --time=05-00:00
#SBATCH -o stdout.%N.%j.out
#SBATCH -e STDERR.%N.%j.err

## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh

for dir in 0_Ir 1_Mn 2_Fe 3_Co 4_Ni; do
    for sub in 1_top_layer 2_hol_layer 3_Ir_top 4_Ir_hol; do
        path="${dir}/${sub}"
        cp $path/* .
	    mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.
3.2.std.x
        cp * $path/
    done
done