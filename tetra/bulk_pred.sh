# #!/bin/bash

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hform base_cfse ICOHPn ionNn mag volume l_bond chg --output period-hform-nn
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output period-hsub-nn

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionNo mag volume l_bond chg --output period-hsub-oo
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNo mag volume l_bond chg --output period-hsub-no

# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_CN
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_OS
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_group
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_outer_e
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag l_bond chg --output wo_volume
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume chg --output wo_l_bond
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn volume l_bond chg --output wo_mag
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond --output wo_chg
# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub ICOHPn ionNn mag volume l_bond chg --output wo_cfse

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionN mag volume l_bond chg --output period-hsub-xx50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionNn mag volume l_bond chg --output period-hsub-xn50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionNo mag volume l_bond chg --output period-hsub-xo50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionNc mag volume l_bond chg --output period-hsub-xc50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionN mag volume l_bond chg --output period-hsub-nx50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionN mag volume l_bond chg --output period-hsub-ox50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionN mag volume l_bond chg --output period-hsub-cx50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output period-hsub-nn50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNo mag volume l_bond chg --output period-hsub-no50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNc mag volume l_bond chg --output period-hsub-nc50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionNn mag volume l_bond chg --output period-hsub-on50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionNo mag volume l_bond chg --output period-hsub-oo50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionNc mag volume l_bond chg --output period-hsub-oc50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionNn mag volume l_bond chg --output period-hsub-cn50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionNo mag volume l_bond chg --output period-hsub-co50
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionNc mag volume l_bond chg --output period-hsub-cc50

python ~/bin/tools/tetra/plot_mae_matrix.py

python ~/bin/tools/tetra/bulk_pred_cfse.py --model rf --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gbr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model xgb --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lgb --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg