# #!/bin/bash

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_OS
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_CN
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn volume l_bond chg --output wo_mag
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag l_bond chg --output wo_volume
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume chg --output wo_l_bond
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond --output wo_chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_group
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_outer_e
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub ICOHPn ionNn mag volume l_bond chg --output wo_cfse

python ~/bin/tools/tetra/bulk_pred_cfse.py --model rf --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gbr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model xgb --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/bulk_pred_cfse.py --model lgb --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg

python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHP ionN mag volume l_bond chg --output hsubxx
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHP ionNn mag volume l_bond chg --output hsubxn
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHP ionNo mag volume l_bond chg --output hsubxo
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHP ionNc mag volume l_bond chg --output hsubxc
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionN mag volume l_bond chg --output hsubnx
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPo ionN mag volume l_bond chg --output hsubox
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPc ionN mag volume l_bond chg --output hsubcx
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output hsubnn
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNo mag volume l_bond chg --output hsubno
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNc mag volume l_bond chg --output hsubnc
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPo ionNn mag volume l_bond chg --output hsubon
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPo ionNo mag volume l_bond chg --output hsuboo
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPo ionNc mag volume l_bond chg --output hsuboc
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPc ionNn mag volume l_bond chg --output hsubcn
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPc ionNo mag volume l_bond chg --output hsubco
python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPc ionNc mag volume l_bond chg --output hsubcc

python ~/bin/tools/tetra/plot_mae_matrix.py

# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hfus base_cfse ICOHPo ionNo mag volume l_bond chg --output Hfus

# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub cfse ICOHPc ionNo mag volume l_bond chg --output cfse

# python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --X OS CN group outer_e Hsub ICOHPc ionNo mag volume l_bond chg --output nocfse
