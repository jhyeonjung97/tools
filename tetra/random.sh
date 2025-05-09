for i in {2..100}; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPo ionNo mag volume l_bond chg --random_state ${i} --output random${i}
done