for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    metal=$(echo "${path[-1]}" | cut -d'_' -f2)

    if [[ -n $(qstat -u x2755a09 | grep "${coord}${row}${numb} ") ]]; then
        if [[ -n $(qstat -u x2755a09 | grep "${coord}${row}${numb}stc") ]]; then
            pwd
            mv opt/* .
            rm -r opt
            qstat -u x2755a09 | grep --color=auto "${coord}${row}${numb}stc"
        fi
    fi
done
