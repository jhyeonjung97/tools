for dir in /pscratch/sd/j/jiuy97/7_V_bulk/*_*_*/5d/04_W/isif2
do
    cd $dir
    python ~/bin/verve/bader.py
done