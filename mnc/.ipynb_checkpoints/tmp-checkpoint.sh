for dir in /pscratch/sd/j/jiuy97/6_MNC/2_OH/*d/*_*/*_*S/*_
do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    # metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    # spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ ! -f ./vasp.out ]]; then
        python ~/bin/tools/mnc/dz.py $dz
    fi
done