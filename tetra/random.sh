for i in {1..100}; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --random_state ${i} --output mendeleev${i} --X OS CN group outer_e Hevap base_cfse ionNo
done

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --output mendeleev --X OS CN group outer_e Hevap base_cfse ionNo
