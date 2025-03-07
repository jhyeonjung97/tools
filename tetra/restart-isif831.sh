i=0
for dir in /pscratch/sd/j/jiuy97/7_V_bulk/6_Octahedral_RS/3d/*_*
do
    cd "$dir" || continue
    
    IFS='/' read -r -a path <<< $dir
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=$(echo "${path[-2]}" | cut -d'_' -f1)
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)

    echo $coord$row$numb
    
    # if [[ -d "isif3" ]]; then
    #     if [[ -f "DONE" ]]; then
    #         continue
    #     elif 
    # elif [[ -d "$dir/isif8/" ]]; then
        
    # elif (( i < 5 )) && [[ -f "$dir/DONE" ]] && [[ -f "$dir/final_opt_bulk8_afm.traj" ]]; then  # Correct integer comparison
    #     pwd
    #     mkdir isif8
    #     #sbatch submit.sh
    #     ((i+=1))  # Correct increment syntax
    # elif (( i < 5 )) && [[ -f "$dir/DONE" ]] && [[ -f "$dir/final_opt_bulk8_afm.traj" ]]; then  # Correct integer comparison
    #     pwd
    #     #sbatch submit.sh
    #     ((i+=1))  # Correct increment syntax
    # fi
done