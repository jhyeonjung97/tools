for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*
do
    cd $dir
    IFS='/' read -r -a path <<< $dir
    vib_dir=/home/hyeonjung/scratch/4_IrFe3/${path[-3]}/vib/${path[-2]}/${path[-1]}
    if [[ -f CONTCAR ]]; then
        python ~/bin/tools/irfe/vib.py
        cp vib.vasp $vib_dir/POSCAR
    fi
    cp ~/bin/tools/irfe/INCAR_vib INCAR
    cp ~/bin/tools/irfe/KPOINTS .
    cp ~/bin/tools/irfe/run_slurm.sh .
done