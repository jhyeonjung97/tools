for i in {1..100}; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output edge${i} \
    --X OS CN group outer_e Hevap base_cfse ICOHPc ionNo mag volume l_bond chg --edge
done

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPc ionNo mag volume l_bond chg --output edge --edge
python ~/bin/tools/tetra/bulk_pred_minus.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPc ionNo mag volume l_bond chg --output minus