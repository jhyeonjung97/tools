for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/*/*_*; do
    if [[ -f $dir/DONE ]] && [[ ! -n $(grep 'PROFILE, used timers') ]]; then
        echo $dir
    fi
done