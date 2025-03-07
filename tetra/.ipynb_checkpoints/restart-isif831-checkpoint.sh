i=0
for dir in /pscratch/sd/j/jiuy97/7_V_bulk/6_Octahedral_RS/3d/*_*
do
    cd "$dir" || continue
    
    IFS='/' read -r -a path <<< $dir
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)

    jobname=${coord}${row}${numb}
    
    if [[ -n $(squeue --me | grep $jobname) ]]; then
        continue
    elif [[ -d "isif3" ]]; then
        if [[ -f "DONE" ]]; then
            continue
        else
            pwd; echo -e "\e[31mCheck this directory!\e[0m"
        fi
    elif [[ -d "isif8" ]]; then
        if [[ -f "DONE" ]]; then
            mkdir isif3; mv * isif3; cd isif3
            mv isif8 $dir; cp CHGCAR WAVECAR restart.json submit.sh $dir; cd $dir
            sed -i -e 's/opt_bulk8/opt_bulk3/' submit.sh
            pwd; echo -e "\e[32mSubmit?\e[0m"; # sbatch submit.sh; ((i+=1))
        else
            pwd; echo -e "\e[31mCheck this directory!\e[0m"
        fi
    else
        if [[ -f "DONE" ]]; then
            mkdir isif8; mv * isif8; cd isif8
            cp CHGCAR WAVECAR restart.json submit.sh $dir; cd $dir
            sed -i -e 's/opt_bulk8/opt_bulk3/' submit.sh
            pwd; echo -e "\e[32mSubmit?\e[0m"; # sbatch submit.sh; ((i+=1))
        else
            pwd; echo -e "\e[31mCheck this directory!\e[0m"
        fi
    fi
done