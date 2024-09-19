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
    sed -i -e "s/jobname/HEO-${metal}${numb}0/" submit.sh
    sbatch submit.sh
    
    mkdir 1_O 2_OH
    cp restart-o.json 1_O/restart.json
    cp restart-oh.json 2_OH/restart.json
    cp ~/bin/tools/heo/submit.sh 1_O/
    cp ~/bin/tools/heo/submit.sh 2_OH/
    sed -i -e "s/jobname/HEO-${metal}${numb}1/" 1_O/submit.sh
    sed -i -e "s/jobname/HEO-${metal}${numb}2/" 2_OH/submit.sh
    cd 1_O; sbatch submit.sh; cd $dir
    cd 2_OH; sbatch submit.sh; cd $dir
done