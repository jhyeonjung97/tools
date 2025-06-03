# #!/bin/bash

# # for coord in 1_Tetrahedral_WZ 2_Tetrahedral_ZB 3_SquarePlanar_TN 4_SquarePlanar_PD 5_SquarePlanar_NB 6_Octahedral_RS 7_Pyramidal_LT 8_Tetrahedral_WZ 9_Tetrahedral_ZB
# # do
# #     for row in 3d 4d 5d
# #     do
# #         for dir in /pscratch/sd/j/jiuy97/8_V_slab/${coord}/${row}/*_*/clean
# #         do
# #             cd ${dir}
# #             if [[ -f 'unmatched' ]] || [[ -f 'DONE' ]]; then
# #                 continue
# #             elif [[ -f 'submit.sh' ]]; then
# #                 sbatch submit.sh
# #             fi
# #         done
# #     done
# # done

for coord in 1_Tetrahedral_WZ 2_Tetrahedral_ZB 3_SquarePlanar_TN 4_SquarePlanar_PD 5_SquarePlanar_NB 6_Octahedral_RS 7_Pyramidal_LT 8_Tetrahedral_WZ 9_Tetrahedral_ZB
do
    for row in 3d 4d 5d
    do
        for ads in clean o oh
        do
            dir="/pscratch/sd/j/jiuy97/8_V_slab/${coord}/${row}/*_*/${ads}"
            cd ${dir}
            IFS='/' read -r -a path <<< $dir
            coordi=$(echo "${path[-3]}" | cut -d'_' -f3)
            numb=$(echo "${path[-1]}" | cut -d'_' -f1)
            jobname=${coordi}${row}${numb}${ads}
            if [[ -f 'unmatched' ]] || [[ -f 'DONE' ]] || [[ -n $(squeue --me | grep $jobname) ]]; then
                continue
            elif [[ -f 'submit.sh' ]]; then
                python ~/bin/get_restart3.py                
            fi
        done
    done
done

# for row in 3d 4d 5d
# do
#     for dir in /pscratch/sd/j/jiuy97/8_V_slab/9_Tetrahedral_ZB/${row}/*_*/clean
#     do
#         cd ${dir}
#         IFS='/' read -r -a path <<< $dir
#         coordi=$(echo "${path[-3]}" | cut -d'_' -f3)
#         numb=$(echo "${path[-1]}" | cut -d'_' -f1)
#         jobname=${coordi}${row}${numb}${ads}
#         if [[ -f 'unmatched' ]] || [[ -f 'DONE' ]]; then
#             continue
#         elif [[ -f 'submit.sh' ]]; then
#             pwd
#         fi
#     done
# done

# for coord in 1_Tetrahedral_WZ 2_Tetrahedral_ZB 5_SquarePlanar_NB 6_Octahedral_RS
# do
#     for row in 3d 4d 5d
#     do
#         for ads in clean o oh
#         do
#             for dir in /pscratch/sd/j/jiuy97/8_V_slab/${coord}/${row}/*_*/${ads}
#             do
#                 cd ${dir}
#                 IFS='/' read -r -a path <<< $dir
#                 coordi=$(echo "${path[-3]}" | cut -d'_' -f3)
#                 numb=$(echo "${path[-1]}" | cut -d'_' -f1)
#                 jobname=${coordi}${row}${numb}${ads}
#                 if [[ -f 'unmatched' ]] || [[ -f 'DONE' ]] || [[ -n $(squeue --me | grep $jobname) ]]; then
#                     continue
#                 elif [[ -f 'submit.sh' ]]; then
#                     python ~/bin/get_restart3
#                 fi
#             done
#         done
#     done
# done


# for row in 3d 4d 5d
# do
#     for dir in /pscratch/sd/j/jiuy97/8_V_slab/9_Tetrahedral_ZB/${row}/*_*/clean
#     do
#         cd ${dir}
#         IFS='/' read -r -a path <<< $dir
#         numb=$(echo "${path[-1]}" | cut -d'_' -f1)
#         jobname=ZZ${row}${numb}
#         if [[ -f 'unmatched' ]] || [[ -f 'DONE' ]] || [[ -n $(squeue --me | grep $jobname) ]]; then
#             continue
#         elif [[ -f 'submit.sh' ]]; then
#             pwd
#         fi
#     done
# done