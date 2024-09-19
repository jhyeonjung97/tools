for dir in /pscratch/sd/j/jiuy97/5_HEO/3_local/*_*/*_
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $dir
    metal=$(echo ${path[-2]} | cut -d'_' -f2)
    numb=$(echo ${path[-1]} | cut -d'_' -f1)
    
    python ~/bin/tools/heo/add-o.py
    python ~/bin/tools/heo/add-oh.py
    
    mkdir o oh
    cp restart-o.json o/restart.json
    cp restart-oh.json oh/restart.json
done

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

for oh_dir in /pscratch/sd/j/jiuy97/6_MNC/2_OH/*_*/*_*S/*_
do
    IFS='/' read -r -a path <<< $oh_dir
    path1=${path[-1]}
    path2=${path[-2]}
    path3=${path[-3]}
    metal=$(echo $path3 | cut -d'_' -f2)
    oh_spin=$(echo $path2 | cut -d'_' -f2)
    dz=$(echo $path1 | cut -d'_' -f1)

    if [[ $metal == 'Mn' ]]; then
        path3='5_Mn'; path4='3d'
    elif [[ $metal == 'Fe' ]]; then
        path3='6_Fe'; path4='3d'
    elif [[ $metal == 'Co' ]]; then
        path3='7_Co'; path4='3d'
    elif [[ $metal == 'Ni' ]]; then
        path3='8_Ni'; path4='3d'
    elif [[ $metal == 'Mo' ]]; then
        path3='4_Mo'; path4='4d'
    elif [[ $metal == 'W' ]]; then
        path3='4_W'; path4='5d'
    fi

    if [[ $metal == 'Mn' || $metal == 'Fe' || $metal == 'Co' || $metal == 'Ni' || $metal == 'Mo' || $metal == 'W' ]]; then
        if [[ $oh_spin == 'LS' ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
        elif [[ $oh_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1" ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1
        elif [[ $oh_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1" ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1
        elif [[ $oh_spin == 'IS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1" ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1
        else
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
        fi
        
        cd $oh_dir; pwd
        cp $dir/restart-oh.json ./restart.json
        cp ~/bin/tools/mnc/submit.sh ./
        sed -i -e "/#SBATCH -J/c\#SBATCH -J OH${metal}${oh_spin}${dz}" submit.sh
        if [[ $oh_spin == 'LS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' submit.sh
        elif [[ $oh_spin == 'IS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' submit.sh
        elif [[ $oh_spin == 'HS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' submit.sh
        fi
        # sbatch submit.sh
    fi
done