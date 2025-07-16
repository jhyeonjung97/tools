for norm in oo no nn on; do
    for feature in Hsub Hfus Hform Hevap; do
        for i in {1..99}; do
            python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output ${feature}-${norm}${i} \
            --X OS CN group outer_e ${feature} base_cfse ICOHPo ionNo mag volume l_bond chg
        done
    done
done