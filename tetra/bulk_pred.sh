# #!/bin/bash

for model in rf lr gpr gbr xgb lgb; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model $model --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
    # python ~/bin/tools/tetra/bulk_pred_edge.py --model $model --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
done

python ~/bin/tools/tetra/bulk_pred_cfse.py --model rf --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lr --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gbr --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model xgb --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lgb --X numb OS CN group outer_e Hevap base_cfse ICOHPmn ionNn mag volume l_bond chg

