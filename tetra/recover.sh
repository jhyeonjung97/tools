# for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*/*_*/; do
#     cd $dir
#     IFS='/' read -r -a path <<< $PWD
#     coord=$(echo "${path[-3]}" | cut -d'_' -f3)
#     row=${path[-2]}
#     numb=$(echo "${path[-1]}" | cut -d'_' -f1)
#     metal=$(echo "${path[-1]}" | cut -d'_' -f2)

#     if [[ -n $(qstat -u x2755a09 | grep "${coord}${row}${numb} ") ]]; then
#         if [[ -n $(qstat -u x2755a09 | grep "${coord}${row}${numb}stc") ]]; then
#             pwd
#             mv opt/* .
#             rm -r opt
#             qstat -u x2755a09 | grep --color=auto "${coord}${row}${numb}stc"
#         fi
#     fi
# done

# for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*/*_*/; do
#     cd $dir
#     IFS='/' read -r -a path <<< $PWD
#     coord=$(echo "${path[-3]}" | cut -d'_' -f3)
#     row=${path[-2]}
#     numb=$(echo "${path[-1]}" | cut -d'_' -f1)
#     metal=$(echo "${path[-1]}" | cut -d'_' -f2)

#     if [[ ! -s lobsterin ]] && [[ -s lobsterout ]]; then
#         cp ~/bin/tools/tetra/lobsterin .
#         sed -i -e "s/X/$metal/g" lobsterin
#         cp ~/bin/tools/tetra/lobster.sh .
#         sed -i -e "s/jobname/${coord}${row}${numb}lob/" lobster.sh
#         pwd
#         qsub lobster.sh
#     elif [[ ! -s icohp.txt ]] && [[ -s lobsterout ]]; then
#         echo -e "\e[36m$PWD\e[0m err"
#     fi
# done

# for dir in /scratch/x2755a09/5_V_bulk/3_SquarePlanar_TN/3d/*_*/; do
#     cd $dir
#     IFS='/' read -r -a path <<< $PWD
#     coord=$(echo "${path[-3]}" | cut -d'_' -f3)
#     row=${path[-2]}
#     numb=$(echo "${path[-1]}" | cut -d'_' -f1)
#     metal=$(echo "${path[-1]}" | cut -d'_' -f2)

#     if [[ -d opt ]] && [[ ! -s icohp.txt ]]; then
#         cp ~/bin/tools/tetra/lobsterin .
#         sed -i -e "s/X/$metal/g" lobsterin
#         cp ~/bin/tools/tetra/static.sh .
#         sed -i -e "s/jobname/${coord}${row}${numb}stc/" static.sh
#         pwd; qsub static.sh
#     fi
# done

for dir in /scratch/x2755a09/5_V_bulk/*_*_*/*d/*_*/; do
    cd $dir
    IFS='/' read -r -a path <<< $PWD
    coord=$(echo "${path[-3]}" | cut -d'_' -f3)
    row=${path[-2]}
    numb=$(echo "${path[-1]}" | cut -d'_' -f1)
    if [[ -n $(qstat -u x2755a09 | grep ${coord}${row}${numb}) ]]; then
        ~/bin/shoulder/rm_mv vasp.out
        ~/bin/shoulder/rm_mv WZ*.e* ZB*.e* TN*.e* PD*.e* NB*.e* RS*.e* LT*.e*
        ~/bin/shoulder/rm_mv WZ*.o* ZB*.o* TN*.o* PD*.o* NB*.o* RS*.o* LT*.o*
    fi
done