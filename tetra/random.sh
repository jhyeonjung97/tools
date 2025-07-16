for i in {1..99}; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output hsub-oo${i} \
    --X OS CN group outer_e Hsub base_cfse ICOHPo ionNo mag volume l_bond chg
done