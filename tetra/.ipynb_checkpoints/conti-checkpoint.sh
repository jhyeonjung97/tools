for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*d/*_*/; do
    cd $dir
    if [[ -d conti_3 ]] && [[ ! -s vasp.out ]]; then
        echo $PWD
        # rm -r conti_3 conti_4 conti_5 conti_6
        # cp conti_2/restart.json .
        # sed -i -e 's/fm.py/fm_ediff.py/' submit.sh
    fi
    # sh ~/bin/verve/conti.sh
    # grep CONTCAR vasp.out
done
