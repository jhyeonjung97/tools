for dir in /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/*d/*_*/
do
    IFS='/' read -r -a path <<< $dir
    path1=${path[-1]}
    path2=${path[-2]}
    path3=${path[-3]}
    numb=$(echo $path1 | cut -d'_' -f1 | cut -c1,2)
    tag=$(echo $path1 | cut -d'_' -f1 | cut -c3)
    metal=$(echo $path1 | cut -d'_' -f2)
    coord=$(echo $path3 | cut -d'_' -f3)

    new_dir="/pscratch/sd/j/jiuy97/3_V_bulk/isif8/${path3}/${path2}/${numb}_${metal}"
    if [[ $tag == 'z' ]] || [[ $tag == 'x' ]] || [[ $tag == 's' ]]; then
        if [[ -s $dir/opt/start.traj ]]; then
            cp $dir/opt/start.traj $new_dir
        else
            cp $dir/opt/POSCAR $new_dir
            rm $new_dir/CONTCAR
        fi
    else
        cp $dir/opt/restart.json $new_dir
    fi
done