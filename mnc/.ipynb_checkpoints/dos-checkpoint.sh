for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*/*_*/*_*/*
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $dir
    row=${path[-4]}
    m=$(echo "${path[-3]}" | cut -d'_' -f1)
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    # if [[ ($metal == 'Mn' || $metal == 'Fe' || $metal == 'Co' || $metal == 'Ni') && ($spin == 'IS' || $spin == 'LS') && ($dz == 0 || $dz == 5) ]]; then
    if [[ $metal == 'Fe' && $spin == 'IS' ]]; then
        if [[ $spin == 'stable' ]]; then
            spin='MS'
        fi
        if [[ $dz == 'relaxed' ]]; then
            dz='r'
        fi
        if [[ ! -s sumo.sh ]]; then
            sed -i 's/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/0.0/g' vasprun.xml
            cp /pscratch/sd/j/jiuy97/6_MNC/scripts/sumo.sh .
            sed -i -e "s/XX/$metal/g" sumo.sh
        fi
        if [[ ! -s sumo_${metal}_dos.dat ]]; then
            sh sumo.sh
        fi
        if [[ ! -s "/pscratch/sd/j/jiuy97/6_MNC/figures/dos/pdos_${row}_${m}${metal}_${spin}${dz}.png" ]]; then
            python ~/bin/tools/mnc/dos.py --file sumo_${metal}_dos.dat --output pdos_${row}_${m}${metal}_${spin}${dz}.png
        fi
    fi
done

for dir in /pscratch/sd/j/jiuy97/6_MNC/*_O*/*_*/*_*/*
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $dir
    ads=$(echo "${path[-4]}" | cut -d'_' -f2)
    m=$(echo "${path[-3]}" | cut -d'_' -f1)
    metal=$(echo "${path[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ $spin == 'stable' ]]; then
        spin='MS'
    fi
    if [[ $dz == 'relaxed' ]]; then
        dz='r'
    fi
    if [[ ! -s sumo.sh ]]; then
        sed -i 's/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/0.0/g' vasprun.xml
        cp /pscratch/sd/j/jiuy97/6_MNC/scripts/sumo.sh .
        sed -i -e "s/XX/$metal/g" sumo.sh
    fi
    if [[ ! -s sumo_${metal}_dos.dat ]]; then
        sh sumo.sh
    fi
    if [[ ! -s "/pscratch/sd/j/jiuy97/6_MNC/figures/dos/pdos_${ads}_${m}${metal}_${spin}${dz}.png" ]]; then
        python ~/bin/tools/mnc/dos.py --file sumo_${metal}_dos.dat --output pdos_${ads}_${m}${metal}_${spin}${dz}.png
    fi
done