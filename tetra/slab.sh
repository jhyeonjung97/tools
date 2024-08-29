qstat -u x2755a09 > ~/mystat.txt

for dir in /scratch/x2755a09/6_V_slab/*_*_*/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)
    
    if [[ ! -n $(grep ${coord}${row}${numb} ~/mystat.txt) ]] && [[ ! -s DONE ]]; then
        echo -e "\e[32m$PWD\e[0m"
    fi
done