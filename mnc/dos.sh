for dir in /pscratch/sd/j/jiuy97/6_MNC/1_O/*_*/*_HS/*_
do
    cd $dir
    IFS='/' read -r -a path <<< $dir
    ads=$(echo "${path[-4]}" | cut -d'_' -f2)
    m=$(echo "${path[-3]}" | cut -d'_' -f1)
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ $spin == 'stable' ]]; then
        spin='MS'
    fi
    if [[ $dz == 0 || $dz == 6 ]] && [[ ! -s "/pscratch/sd/j/jiuy97/6_MNC/figures/pdos/pdos_${ads}_${m}${metal}_${spin}${dz}.png" ]]; then
        sed -i 's/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/0.0/g' vasprun.xml
        cp /pscratch/sd/j/jiuy97/6_MNC/scripts/sumo.sh .
        sed -i -e "s/XX/$metal/g" sumo.sh
        sh sumo.sh
        python ~/bin/tools/mnc/dos.py --file sumo_${metal}_dos.dat --output pdos_${ads}_${m}${metal}_${spin}${dz}.png
    fi
done