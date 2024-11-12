for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/*_*/*_*S/
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< "$dir"
    metal=$(echo "${path[-2]}" | cut -d'_' -f2)
    spin=$(echo "${path[-1]}" | cut -d'_' -f2)
    for ads in o oh
    do
        if [[ -f $ads/DONE ]] && [[ ! -d $ads/icohp ]] && [[ ! $spin == 'IS' ]]; then
            cd $ads
            mkdir -p icohp
            cp WAVECAR restart.json icohp/
            cp ~/bin/tools/mnc/lobsterin icohp/
            sed -i -e "s/X/${metal}/g" icohp/lobsterin
            cp ~/bin/tools/mnc/static.sh icohp/submit.sh
            sed -i -e "s/jobname/${metal}${ads}COHP/g" icohp/submit.sh
            cd icohp/
            sbatch submit.sh
            cd $dir
        else
            continue
        fi
    done
done