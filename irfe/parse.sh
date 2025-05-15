# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*
# do
#     cd $dir; pwd
#     if [[ -f OUTCAR ]] && [[ ! -f final_with_calculator.json ]]; then
#         ase convert -n -1 OUTCAR final_with_calculator.json
#     fi
#     if [[ -f AECCAR0 ]] && [[ ! -f atoms_bader_charge.json ]]; then
#         python3 ~/bin/verve/bader.py
#     fi
# done

# for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/vib/*_*/*_*_*
# do
#     cd $dir; pwd
#     IFS='/' read -r -a path <<< $dir
#     if [[ ${path[-1]} == "5_M_top" ]] || [[ ${path[-1]} == "6_M_hol" ]]; then
#         continue
#     fi
#     if [[ -f OUTCAR ]] && [[ ! -f vib.txt ]]; then
#         vaspkit -task 501 > vib.txt
#     fi
# done

for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*/
do
    cd $dir; pwd
    if [[ -f DONE ]]; then
        ase convert -f -n -1 OUTCAR final_with_calculator.json
    else
        rm final_with_calculator.json
    fi
    if [[ -f AECCAR0 ]]; then
        python3 ~/bin/verve/bader.py
    fi
    if [[ -f $dir/vib/OUTCAR ]]; then
        cd $dir/vib
        vaspkit -task 501 > vib.txt
        cp vib.txt $dir/
    fi
done