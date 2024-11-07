# for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*S/*_
# do
#     cd $dir; pwd
#     python ~/bin/tools/mnc/add-o.py
#     python ~/bin/tools/mnc/add-oh.py
#     python ~/bin/tools/mnc/add-ooh.py
#     python ~/bin/tools/mnc/add-co.py
#     python ~/bin/tools/mnc/add-h.py
# done

# for o_dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/*_*/*_*S/*_
# do
#     IFS='/' read -r -a path <<< $o_dir
#     path1=${path[-1]}
#     path2=${path[-2]}
#     path3=${path[-3]}
#     metal=$(echo $path3 | cut -d'_' -f2)
#     o_spin=$(echo $path2 | cut -d'_' -f2)
#     dz=$(echo $path1 | cut -d'_' -f1)

#     if [[ $metal == 'Mn' ]]; then
#         path3='5_Mn'; path4='3d'
#     elif [[ $metal == 'Fe' ]]; then
#         path3='6_Fe'; path4='3d'
#     elif [[ $metal == 'Co' ]]; then
#         path3='7_Co'; path4='3d'
#     elif [[ $metal == 'Ni' ]]; then
#         path3='8_Ni'; path4='3d'
#     elif [[ $metal == 'Mo' ]]; then
#         path3='4_Mo'; path4='4d'
#     elif [[ $metal == 'W' ]]; then
#         path3='4_W'; path4='5d'
#     fi

#     if [[ $metal == 'Mn' || $metal == 'Fe' || $metal == 'Co' || $metal == 'Ni' || $metal == 'Mo' || $metal == 'W' ]]; then
#         if [[ $o_spin == 'LS' ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
#         elif [[ $o_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1
#         elif [[ $o_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1
#         elif [[ $o_spin == 'IS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1
#         else
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
#         fi
        
#         cd $o_dir; pwd
#         cp $dir/restart-o.json ./restart.json
#         cp $dir/WAVECAR ./
#         cp ~/bin/tools/mnc/submit.sh ./
#         sed -i -e "/#SBATCH -J/c\#SBATCH -J O${metal}${o_spin}${dz}" submit.sh
#         if [[ $o_spin == 'LS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' submit.sh
#         elif [[ $o_spin == 'IS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' submit.sh
#         elif [[ $o_spin == 'HS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' submit.sh
#         fi
#         # sbatch submit.sh
#     fi
# done

# for oh_dir in /pscratch/sd/j/jiuy97/6_MNC/2_OH/*_*/*_*S/*_
# do
#     IFS='/' read -r -a path <<< $oh_dir
#     path1=${path[-1]}
#     path2=${path[-2]}
#     path3=${path[-3]}
#     metal=$(echo $path3 | cut -d'_' -f2)
#     oh_spin=$(echo $path2 | cut -d'_' -f2)
#     dz=$(echo $path1 | cut -d'_' -f1)

#     if [[ $metal == 'Mn' ]]; then
#         path3='5_Mn'; path4='3d'
#     elif [[ $metal == 'Fe' ]]; then
#         path3='6_Fe'; path4='3d'
#     elif [[ $metal == 'Co' ]]; then
#         path3='7_Co'; path4='3d'
#     elif [[ $metal == 'Ni' ]]; then
#         path3='8_Ni'; path4='3d'
#     elif [[ $metal == 'Mo' ]]; then
#         path3='4_Mo'; path4='4d'
#     elif [[ $metal == 'W' ]]; then
#         path3='4_W'; path4='5d'
#     fi

#     if [[ $metal == 'Mn' || $metal == 'Fe' || $metal == 'Co' || $metal == 'Ni' || $metal == 'Mo' || $metal == 'W' ]]; then
#         if [[ $oh_spin == 'LS' ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
#         elif [[ $oh_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1
#         elif [[ $oh_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1
#         elif [[ $oh_spin == 'IS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1
#         else
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
#         fi
        
#         cd $oh_dir; pwd
#         cp $dir/restart-oh.json ./restart.json
#         cp ~/bin/tools/mnc/submit.sh ./
#         sed -i -e "/#SBATCH -J/c\#SBATCH -J OH${metal}${oh_spin}${dz}" submit.sh
#         if [[ $oh_spin == 'LS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' submit.sh
#         elif [[ $oh_spin == 'IS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' submit.sh
#         elif [[ $oh_spin == 'HS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' submit.sh
#         fi
#         # sbatch submit.sh
#     fi
# done

# for ooh_dir in /pscratch/sd/j/jiuy97/6_MNC/3_OOH/*_*/*_*S/*_
# do
#     IFS='/' read -r -a path <<< $ooh_dir
#     path1=${path[-1]}
#     path2=${path[-2]}
#     path3=${path[-3]}
#     metal=$(echo $path3 | cut -d'_' -f2)
#     ooh_spin=$(echo $path2 | cut -d'_' -f2)
#     dz=$(echo $path1 | cut -d'_' -f1)

#     if [[ $metal == 'Mn' ]]; then
#         path3='5_Mn'; path4='3d'
#     elif [[ $metal == 'Fe' ]]; then
#         path3='6_Fe'; path4='3d'
#     elif [[ $metal == 'Co' ]]; then
#         path3='7_Co'; path4='3d'
#     elif [[ $metal == 'Ni' ]]; then
#         path3='8_Ni'; path4='3d'
#     elif [[ $metal == 'Mo' ]]; then
#         path3='4_Mo'; path4='4d'
#     elif [[ $metal == 'W' ]]; then
#         path3='4_W'; path4='5d'
#     fi

#     if [[ $metal == 'Mn' || $metal == 'Fe' || $metal == 'Co' || $metal == 'Ni' || $metal == 'Mo' || $metal == 'W' ]]; then
#         if [[ $ooh_spin == 'LS' ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
#         elif [[ $ooh_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1
#         elif [[ $ooh_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1
#         elif [[ $ooh_spin == 'IS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1" ]]; then
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1
#         else
#             dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
#         fi
        
#         cd $ooh_dir; pwd
#         cp $dir/restart-ooh.json ./restart.json
#         cp ~/bin/tools/mnc/submit.sh ./
#         sed -i -e "/#SBATCH -J/c\#SBATCH -J OOH${metal}${ooh_spin}${dz}" submit.sh
#         if [[ $ooh_spin == 'LS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' submit.sh
#         elif [[ $ooh_spin == 'IS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' submit.sh
#         elif [[ $ooh_spin == 'HS' ]]; then
#             sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' submit.sh
#         fi
#         # sbatch submit.sh
#     fi
# done

for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/*_*/*_*S
do
    cd "$dir" || continue
    pwd
    IFS='/' read -r -a path <<< "$dir"
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    [[ $spin == 'IS' ]] && continue
    
    cd relaxed || continue
    python ~/bin/tools/mnc/add-o.py
    python ~/bin/tools/mnc/add-oh.py
    # python ~/bin/tools/mnc/add-ooh.py
    cd "$dir" || continue
    
    for ads in o oh #ooh
    do
        mkdir -p "$ads"
        cd "$ads" || exit
        mv ../relaxed/restart-${ads}.json restart.json
        ase convert restart.json POSCAR
        ase convert POSCAR start.traj
        rm restart.json
        
        cp ~/bin/tools/mnc/submit.sh ./
        sed -i -e "/#SBATCH -J/c\#SBATCH -J ${ads}${metal}${spin}r" submit.sh
        if [[ $spin == 'LS' ]]; then
            sed -i -e "s/mnc-sol.py/mnc-sol-ls.py/" submit.sh
        elif [[ $spin == 'HS' ]]; then
            sed -i -e "s/mnc-sol.py/mnc-sol-hs.py/" submit.sh
        fi
        # sbatch submit.sh
        cd "$dir"
    done

    # mkdir -p ooh
    # cd ooh || exit
    # mv ../relaxed/restart-ooh.json restart.json
    # cp ~/bin/tools/mnc/submit.sh ./
    # sed -i -e "/#SBATCH -J/c\#SBATCH -J ooh${metal}MSr" submit.sh
    # sbatch submit.sh
done
