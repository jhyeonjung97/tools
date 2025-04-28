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

# for dir in 5_IrMn 6_IrFe 7_IrCo 8_IrNi
# do
#     cd /home/hyeonjung/scratch/4_IrFe3/slab/$dir
#     cp h-top-atom1.vasp /home/hyeonjung/scratch/4_IrFe3/1_H/$dir/3_atom_top1
#     cp h-hol-atom1.vasp /home/hyeonjung/scratch/4_IrFe3/1_H/$dir/5_atom_hol1
#     cp o-top-atom1.vasp /home/hyeonjung/scratch/4_IrFe3/3_O/$dir/3_atom_top1
#     cp o-hol-atom1.vasp /home/hyeonjung/scratch/4_IrFe3/3_O/$dir/5_atom_hol1

#     for subdir in /home/hyeonjung/scratch/4_IrFe3/*_*/$dir/*_*_*1
#     do
#         cd $subdir
#         if [[ ! -f OUTCAR ]]; then
#             cp ~/bin/tools/irfe/INCAR .
#             cp ~/bin/tools/irfe/KPOINTS .
#             mv *.vasp POSCAR
#             vaspkit -task 107
#             mv POSCAR_REV POSCAR
#             rm POTCAR
#             vaspkit -task 103
#             python3 ~/bin/orange/magmom.py
#         fi
#     done
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/5_OXRb/*_*_*/*_*_*
# do
#     cd $dir
#     if [[ -f *.vasp ]]; then
#         pwd
#             cp ~/bin/tools/irfe/INCAR .
#             cp ~/bin/tools/irfe/KPOINTS .
#             mv *.vasp POSCAR
#             vaspkit -task 107
#             mv POSCAR_REV POSCAR
#             rm POTCAR
#             vaspkit -task 103
#             python3 ~/bin/orange/magmom.py
#     fi
# done