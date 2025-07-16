for i in {1..99}; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output random${i} --X OS CN group outer_e Hform base_cfse ICOHPo ionNo mag volume l_bond chg
done