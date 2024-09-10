for dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/*d/*_*/*_*S/*_/
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $PWD
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ $dz == '2' ]] || [[ $dz == '4' ]] || [[ $dz == '6' ]]; then
        sbatch submit.sh
    fi
done