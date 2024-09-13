# for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*S/*_
# do
#     cd $dir; pwd
#     python ~/bin/tools/mnc/add-o.py
#     python ~/bin/tools/mnc/add-oh.py
#     python ~/bin/tools/mnc/add-co.py
#     python ~/bin/tools/mnc/add-h.py
# done

for o_dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/*_*/*_*S/*_
do
    IFS='/' read -r -a path <<< $o_dir
    path1=${path[-1]}
    path2=${path[-2]}
    path3=${path[-3]}
    metal=$(echo $path3 | cut -d'_' -f2)
    o_spin=$(echo $path2 | cut -d'_' -f2)
    dz=$(echo $path1 | cut -d'_' -f1)

    if [[ $metal == 'Mn' || $metal == 'Fe' || $metal == 'Co' || $metal == 'Ni' ]]; then
        path4='3d'
    elif [[ $metal == 'Mo' ]]; then
        path4='4d'
    elif [[ $metal == 'W' ]]; then
        path4='5d'
    fi

    if [[ $metal == 'Mn' || $metal == 'Fe' || $metal == 'Co' || $metal == 'Ni' || $metal == 'Mo' || $metal == 'W' ]]; then
        if [[ $o_spin == 'LS' ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
        elif [[ $o_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1" ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_HS/$path1
        elif [[ $o_spin == 'HS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1" ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/3_HS/$path1
        elif [[ $o_spin == 'IS' ]] && [[ -d "/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1" ]]; then
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/2_IS/$path1
        else
            dir=/pscratch/sd/j/jiuy97/6_MNC/0_clean/$path4/$path3/1_LS/$path1
        fi
        
        cd $o_dir; pwd
        cp $dir/restart_o.json ./restart.json
        cp $dir/WAVECAR ./
        cp ~/bin/tools/mnc/submit.sh ./
        sed -i -e "/#SBATCH -J/c\#SBATCH -J O${metal}${o_spin}${dz}" submit.sh
        if [[ $o_spin == 'LS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-ls-nupdown.py/' ${dz}_/submit.sh
        elif [[ $o_spin == 'IS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-is-nupdown.py/' ${dz}_/submit.sh
        elif [[ $o_spin == 'HS' ]]; then
            sed -i 's/mnc-sol.py/mnc-sol-hs-nupdown.py/' ${dz}_/submit.sh
        fi
        # sbatch submit.sh
    fi
done