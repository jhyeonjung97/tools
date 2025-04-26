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

# for dir in 6_IrFe 7_IrCo 8_IrNi # 5_IrMn
# do
#     cd /home/hyeonjung/scratch/4_IrFe3/slab/$dir
#     cp h-*.vasp /home/hyeonjung/scratch/4_IrFe3/1_H/$dir
#     cp o-*.vasp /home/hyeonjung/scratch/4_IrFe3/3_O/$dir
#     cp oh-*.vasp /home/hyeonjung/scratch/4_IrFe3/2_OH/$dir

#     cd /home/hyeonjung/scratch/4_IrFe3/1_H/$dir
#     mv h-top-layer.vasp 1_layer_top
#     mv h-hol-layer.vasp 2_layer_hol
#     mv h-top-atom.vasp 3_atom_top
#     mv h-hol-atom.vasp 4_atom_hol

#     cd /home/hyeonjung/scratch/4_IrFe3/2_OH/$dir
#     mv oh-top-layer.vasp 1_layer_top
#     mv oh-brg-layer.vasp 2_layer_brg
#     mv oh-top-atom1.vasp 3_atom_top1
#     mv oh-top-atom2.vasp 4_atom_top2
#     mv oh-brg-atom1.vasp 5_atom_brg1
#     mv oh-brg-atom2.vasp 6_atom_brg2

#     cd /home/hyeonjung/scratch/4_IrFe3/3_O/$dir
#     mv o-top-layer.vasp 1_layer_top
#     mv o-hol-layer.vasp 2_layer_hol
#     mv o-top-atom.vasp 3_atom_top
#     mv o-hol-atom.vasp 4_atom_hol

#     for subdir in /home/hyeonjung/scratch/4_IrFe3/*_*/$dir/*_*_*
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

for dir in /home/hyeonjung/scratch/4_IrFe3/5_OXRb/*_*_*/*_*_*
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