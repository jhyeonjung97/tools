for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/*/*_*; do
    if [[ ! -f $dir/start.traj ]]; then
        echo $dir; rm $dir/*
    fi
done