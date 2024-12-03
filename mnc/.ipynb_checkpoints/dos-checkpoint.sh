for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/*_*/*_
do
    cd $dir
    IFS='/' read -r -a path <<< $dir
    row=${path[-4]}
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ $metal == 'Ti' ]] && [[ $spin == 'HS' ]] && [[ $dz == '2' ]]; then
        sed -i 's/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/0.0/g' vasprun.xml
        cp /pscratch/sd/j/jiuy97/6_MNC/scripts/sumo.sh .
        sed -i -e "s/XX/$metal/g" sumo.sh
        sh sumo.sh
        python ~/bin/tools/mnc/dos.py --file sumo_${metal}_dos.dat --output pdos_${metal}${spin}${dz}.png
    fi
done