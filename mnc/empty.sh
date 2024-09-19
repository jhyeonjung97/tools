for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*
do
    cd $dir; pwd
    mkdir empty; cd empty
    python ~/bin/tools/mnc/empty.py
    cp ~/bin/tools/mnc/submit.sh .
    sed -i -e "s/jobname/${metal}e/"
    
done
