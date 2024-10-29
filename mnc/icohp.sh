for dir in /pscratch/sd/j/jiuy97/6_MNC/0_clean/*d/*_*/most_stable
do
    cd $dir; pwd
    for ads in o oh
    do
        if [[ -f $ads/DONE ]]; then
            cd $ads
            
        else
            exit
        fi
    done
done