for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*/*_*/*_*S;
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)

    path3=${path[-3]}
    path2=${path[-2]}
    path1=${path[-1]}
    ls -d /pscratch/sd/j/jiuy97/6_MNC/kisti/3_MNC/0_clean/$path3/$path2/$path1
done