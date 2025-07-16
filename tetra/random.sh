for norm1 in o n; do
    for norm2 in o n; do
        for feature in Hsub Hfus Hform Hevap; do
            for i in {1..99}; do
                python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output ${feature}-${norm1}${norm2}${i} \
                --X OS CN group outer_e ${feature} base_cfse ICOHP${norm1} ion${norm2} mag volume l_bond chg
            done
        done
    done
done