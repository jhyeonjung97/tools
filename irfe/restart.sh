# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/5_M_top
# do
#     cd $dir
#     find . ! -name '*top*.vasp' -delete
#     cp *.vasp POSCAR
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/6_M_hol
# do
#     cd $dir
#     find . ! -name '*hol*.vasp' -delete
#     cp *.vasp POSCAR
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_M_*
# do
#     cd $dir
#     cp ~/bin/tools/irfe/INCAR .
#     cp ~/bin/tools/irfe/KPOINTS .
#     cp *.vasp POSCAR
#     vaspkit -task 107
#     mv POSCAR_REV POSCAR
#     rm POTCAR
#     vaspkit -task 103
#     python3 ~/bin/orange/magmom.py
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*
# do
#     cd $dir
#     cp ~/bin/tools/irfe/run_slurm.sh .
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/2_OH/*_*/*_*_*
# do
#     cd $dir
#     if [[ ! -f OUTCAR ]]; then
#         cp ~/bin/tools/irfe/INCAR .
#         cp ~/bin/tools/irfe/KPOINTS .
#         cp *.vasp POSCAR
#         vaspkit -task 107
#         mv POSCAR_REV POSCAR
#         rm POTCAR
#         vaspkit -task 103
#         python3 ~/bin/orange/magmom.py
#     fi
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/*_Ir*/*_*_*/*_*_*
# do
#     cd $dir
#     if [[ ! -f OUTCAR ]]; then
#         cp ~/bin/tools/irfe/INCAR .
#         cp ~/bin/tools/irfe/KPOINTS .
#         mv *.vasp POSCAR
#         vaspkit -task 107
#         mv POSCAR_REV POSCAR
#         rm POTCAR
#         vaspkit -task 103
#         python3 ~/bin/orange/magmom.py
#     fi
# done

for dir in /home/hyeonjung/scratch/4_IrFe3/2_OH/*_*/*_*_brg
do
    cd $dir
    if [[ ! -f OUTCAR ]]; then
        cp ~/bin/tools/irfe/INCAR .
        cp ~/bin/tools/irfe/KPOINTS .
        mv *.vasp POSCAR
        vaspkit -task 107
        mv POSCAR_REV POSCAR
        rm POTCAR
        vaspkit -task 103
        python3 ~/bin/orange/magmom.py
    fi
done