for norm1 in o n; do
    for norm2 in o n; do
        for feature in Hsub Hform; do
            for i in {1..99}; do
                if [[ -f ~/Desktop/7_V_bulk/figures/form_pred_cfse_gpr00_all_${feature}${norm1}${norm2}${i}.log ]]; then
                    continue
                fi
                python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output ${feature}${norm1}${norm2}${i} \
                --X OS CN group outer_e ${feature} base_cfse ICOHP${norm1} ion${norm2} mag volume l_bond chg
            done
        done
    done
done