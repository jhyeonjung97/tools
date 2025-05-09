for i in {2..100}; do
    for norm_a in c o n; do
        for norm_b in c o n; do
            echo "\033[1;31mrandom_state: ${i}, norm_a: ${norm_a}, norm_b: ${norm_b}\033[0m"
            name_a=$(echo ${norm_a} | tr '[:lower:]' '[:upper:]')
            name_b=$(echo ${norm_b} | tr '[:lower:]' '[:upper:]')
            python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHP${norm_a} ionN${norm_b} mag volume l_bond chg --random_state ${i} --output random${name_a}${name_b}${i}
            python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap ICOHP${norm_a} ionN${norm_b} mag volume l_bond chg --random_state ${i} --output random${name_a}${name_b}${i}_nocfse
        done
    done
done