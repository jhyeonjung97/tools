#!/bin/bash

# bulkox bulk_e bulk_ie bulk_ieo bulk_ien bulk_cfse
# surfox surf_e surf_ie surf_ieo surf_ien surf_cfse
# ICOHP ICOHPo ICOHPn Hform Hsub Hevap Hfus
# CN group period chg mag volume l_bond cfse

python ~/bin/tools/tetra/lr-slab.py --Y o_energy --output compact \
--X bulkox bulk_e bulk_ieo surfox surf_e surf_ieo CN group period ICOHPo chg mag cfse volume l_bond

python ~/bin/tools/tetra/slab_pred.py --Y oh \
--X bulkox bulk_e bulk_ieo bulk_cfse surfox surf_e surf_ieo surf_cfse CN group period ICOHPo chg mag cfse volume l_bond

python ~/bin/tools/tetra/slab_pred.py --Y oh --output comer \
--X bulkox bulk_e bulk_ieo surfox surf_e surf_ieo group ICOHPo

python ~/bin/tools/tetra/slab_pred_nn.py --Y oh --output comer \
--X bulkox bulk_e bulk_ieo surfox surf_e surf_ieo group ICOHPo

CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg period --output w_period
python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output wo_CN
python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group outer_e Hsub ICOHPn ionNn mag volume l_bond chg --output wo_cfse
python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn volume l_bond chg --output wo_mag
python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag l_bond chg --output wo_volume
python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume chg --output wo_l_bond
python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond --output wo_chg

# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionN mag volume l_bond chg --output xx50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionNn mag volume l_bond chg --output xn50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionNo mag volume l_bond chg --output xo50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHP ionNc mag volume l_bond chg --output xc50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionN mag volume l_bond chg --output nx50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionN mag volume l_bond chg --output ox50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionN mag volume l_bond chg --output cx50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg --output nn50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNo mag volume l_bond chg --output no50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPn ionNc mag volume l_bond chg --output nc50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionNn mag volume l_bond chg --output on50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionNo mag volume l_bond chg --output oo50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPo ionNc mag volume l_bond chg --output oc50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionNn mag volume l_bond chg --output cn50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionNo mag volume l_bond chg --output co50
# python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group period outer_e Hsub base_cfse ICOHPc ionNc mag volume l_bond chg --output cc50
# python ~/bin/tools/tetra/plot_mae_matrix.py

python ~/bin/tools/tetra/slab_pred.py --model rf --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/slab_pred.py --model lr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/slab_pred.py --model gpr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/slab_pred.py --model gbr --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/slab_pred.py --model xgb --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg
python ~/bin/tools/tetra/slab_pred.py --model lgb --X OS CN group outer_e Hsub base_cfse ICOHPn ionNn mag volume l_bond chg