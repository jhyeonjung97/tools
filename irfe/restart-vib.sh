# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*
# do
#     cd $dir
#     IFS='/' read -r -a path <<< $dir
#     vib_dir=/home/hyeonjung/scratch/4_IrFe3/${path[-3]}/vib/${path[-2]}/${path[-1]}
#     if [[ -f CONTCAR ]]; then
#         python ~/bin/tools/irfe/vib.py
#         cp vib.vasp $vib_dir/POSCAR
#     fi
#     cd $vib_dir
#     cp ~/bin/tools/irfe/INCAR_vib INCAR
#     cp ~/bin/tools/irfe/KPOINTS .
#     vaspkit -task 107
#     mv POSCAR_REV POSCAR
#     vaspkit -task 103
#     python3 ~/bin/orange/magmom.py
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/6_M_hol
# do
#     cd $dir; pwd
#     IFS='/' read -r -a path <<< $dir
#     vib_dir=/home/hyeonjung/scratch/4_IrFe3/${path[-3]}/vib/${path[-2]}/${path[-1]}
#     if [[ -f CONTCAR ]]; then
#         python ~/bin/tools/irfe/vib.py
#         cp vib.vasp $vib_dir/POSCAR
#     fi
#     cd $vib_dir
#     cp ~/bin/tools/irfe/INCAR_vib INCAR
#     cp ~/bin/tools/irfe/KPOINTS .
#     vaspkit -task 107
#     mv POSCAR_REV POSCAR
#     vaspkit -task 103
#     python3 ~/bin/orange/magmom.py
# done

for dir in /home/hyeonjung/scratch/4_IrFe3/sgr/*_*_*/*_*
do
    cd $dir
    if [[ -f CONTCAR ]]; then
        python ~/bin/tools/irfe/vib.py
        mkdir -p vib
        cp vib.vasp vib/POSCAR
    fi
    cd vib
    cp ~/bin/tools/irfe/INCAR_vib INCAR
    cp ~/bin/tools/irfe/KPOINTS .
    vaspkit -task 107
    mv POSCAR_REV POSCAR
    vaspkit -task 103
    python3 ~/bin/orange/magmom.py
done

# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*
# do
#     cd $dir
#     if [[ -d vib ]] && [[ -f unmatched ]]; then
#         rm -r vib
#     fi
# done