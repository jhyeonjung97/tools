for dir in /pscratch/sd/j/jiuy97/6_MNC/icohp/*_*/*_/
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< "$dir"
    ads=$(echo "${path[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path[-1]}" | cut -d'_' -f1)
    cp ~/bin/tools/mnc/lobsterin2 ./lobsterin    
    cp ~/bin/tools/mnc/static.sh ./submit.sh
    sed -i -e "s/jobname/ICOHP_${ads}${dz}/" icohp/submit.sh
    sbatch submit.sh
done