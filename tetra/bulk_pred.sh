# #!/bin/bash

for model in rf lr gpr gbr xgb lgb; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model $model --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
    python ~/bin/tools/tetra/bulk_pred_edge.py --model $model --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
done

python ~/bin/tools/tetra/bulk_pred_cfse.py --model rf --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lr --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gbr --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model xgb --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lgb --X OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPn ionNn mag volume l_bond chg --output normnn
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPn ionNo mag volume l_bond chg --output normno
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPn ionNc mag volume l_bond chg --output normnc
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPo ionNn mag volume l_bond chg --output normon
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPo ionNo mag volume l_bond chg --output normoo
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPo ionNc mag volume l_bond chg --output normoc
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPc ionNn mag volume l_bond chg --output normcn
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPc ionNo mag volume l_bond chg --output normco
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPc ionNc mag volume l_bond chg --output normcc
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHP ionN mag volume l_bond chg --output normxx


python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hevap base_cfse ICOHPo ionNo mag volume l_bond chg