#!/bin/bash

for model in rf lr gpr gbr xgb lgb; do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model $model --X numb OS chg l_bond ICOHPmn ion+1n ICOOPmn Natom d_electrons base_cfse exchange_stabilization CN mag madelung
    python ~/bin/tools/tetra/bulk_pred_edge.py --model $model --X numb OS chg l_bond ICOHPmn ion+1n ICOOPmn Natom d_electrons base_cfse exchange_stabilization CN mag madelung
done

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X numb OS chg l_bond ICOHPmn ion+1n ICOOPmn Natom d_electrons base_cfse exchange_stabilization CN mag madelung
python ~/bin/tools/tetra/bulk_pred_edge.py --model gpr --X numb OS chg l_bond ICOHPmn ion+1n ICOOPmn Natom d_electrons base_cfse exchange_stabilization CN mag madelung
