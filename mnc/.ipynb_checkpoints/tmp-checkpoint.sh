for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*/*_*/*_*S;
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    path3=${path[-3]}
    path2=${path[-2]}
    path1=${path[-1]}
    mv $dir/nupdown/* /pscratch/sd/j/jiuy97/6_MNC/kisti/3_MNC/0_clean/$path3/$path2/$path1/0_
done