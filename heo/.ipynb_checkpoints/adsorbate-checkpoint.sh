for dir in /pscratch/sd/j/jiuy97/5_HEO/3_local/*_*/*_
do
    cd $dir; pwd
    IFS='/' read -r -a path <<< $dir
    metal=$(echo ${path[-2]} | cut -d'_' -f2)
    numb=$(echo ${path[-1]} | cut -d'_' -f1)
    
    python ~/bin/get_restart3
    python ~/bin/tools/heo/add-o.py
    python ~/bin/tools/heo/add-oh.py
    
    cp ~/bin/tools/heo/submit.sh .
    sed -i -e "s/jobname/HEO-${metal}${numb}/" submit.sh
    # sbatch submit.sh
    
    mkdir o oh
    cp restart-o.json o/restart.json
    cp restart-oh.json oh/restart.json
    cp submit.sh o/
    cp submit.sh oh/
    sed -i -e "s/jobname/HEO-${metal}${numb}/" submit.sh
    sed -i -e "s/jobname/HEO-${metal}${numb}/" submit.sh
    # cd o/; sbatch submit.sh; cd $dir
    # cd oh/; sbatch submit.sh; cd $dir
done