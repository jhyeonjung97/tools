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


# for dir in /home/hyeonjung/scratch/4_IrFe3/2_OH/5_IrMn/*_*_*
# do
#     cd $dir
#     cp CONTCAR POSCAR
#     cp ~/bin/tools/irfe/INCAR .
#     cp ~/bin/tools/irfe/KPOINTS .
#     vaspkit -task 107
#     mv POSCAR_REV POSCAR
#     rm POTCAR
#     vaspkit -task 103
#     python3 ~/bin/orange/magmom.py
# done
# cd ../..
# ls

for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*
do
    IFS='/' read -r -a path <<< $dir
    site=${path[-1]}
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    ads=$(echo "${path[-3]}" | cut -d'_' -f2)

    if [[ $ads == 'H' ]]; then
        if [[ $metal == 'Fe' ]] || [[ $metal == 'Co' ]] || [[ $metal == 'Ni' ]] || [[ $metal == 'Mn' ]]; then
            if [[ $site == '1_layer_top' ]] || [[ $site == '5_M_top' ]]; then
                continue
            fi
        fi
        if [[ $metal == 'IrFe' ]] || [[ $metal == 'IrCo' ]] || [[ $metal == 'IrNi' ]] || [[ $metal == 'IrMn' ]]; then
            if [[ $site == '2_layer_hol' ]] || [[ $site == '6_atom_hol2' ]]; then
                continue
            fi
        fi
    elif [[ $ads == 'OH' ]]; then
        if [[ $metal == 'Fe' ]] || [[ $metal == 'Co' ]] || [[ $metal == 'Ni' ]] || [[ $metal == 'Mn' ]]; then
            if [[ $site == '1_layer_top' ]] || [[ $site == '5_M_top' ]]; then
                continue
            fi
        fi
        if [[ $metal == 'IrFe' ]] || [[ $metal == 'IrCo' ]] || [[ $metal == 'IrNi' ]] || [[ $metal == 'IrMn' ]]; then
            if [[ $site == '2_layer_brg' ]]; then
                continue
            fi
        fi
    elif [[ $ads == 'O' ]]; then
        if [[ $metal == 'Fe' ]] || [[ $metal == 'Co' ]] || [[ $metal == 'Ni' ]] || [[ $metal == 'Mn' ]]; then
            if [[ $site == '5_M_top' ]]; then
                continue
            fi
        fi
        if [[ $metal == 'IrFe' ]] || [[ $metal == 'IrCo' ]] || [[ $metal == 'IrNi' ]] || [[ $metal == 'IrMn' ]]; then
            if [[ $site == '2_layer_hol' ]]; then
                continue
            fi
        fi
    elif [[ $ads == 'R1' ]] || [[ $ads == 'R2' ]]; then
        site=${path[-2]}
        if [[ $site == '6_IrFe_hol1' ]] || [[ $site == '7_IrFe_hol2' ]] || [[ $site == '8_IrFe_hol3' ]]; then
            continue
        fi
    else
        continue
    fi

    echo $ads $metal $site
done