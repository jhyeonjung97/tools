for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*d/*_*/; do
    cd $dir
    if [[ -n $(grep CONTCAR vasp.out) ]]; then
        if [[ -d conti_2 ]]; then
            ls
        fi
        # grep CONTCAR vasp.out
    fi
done
