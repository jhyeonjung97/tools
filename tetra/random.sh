for i in {1..100}; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output edge${i} \
    --X OS CN group outer_e Hevap base_cfse ICOHPc ionNo mag volume l_bond chg --filter_n_electrons
done