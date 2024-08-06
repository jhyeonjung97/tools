for dir in /scratch/x2755a09/5_V_bulk/*_*_*/fm/*_*/; do
    cd $dir
    if [[ -n $(grep CONTCAR vasp.out) ]]; then
        echo $PWD
        if [[ -d conti_2 ]]; then
            rm -r conti_3 conti_4 conti_5 conti_6
            cp conti_2/* .
            sed -i -e 's/fm.py/fm_ediff.py' submit.sh
        fi
        sh ~/bin/verve/conti.sh
        # grep CONTCAR vasp.out
    fi
done
