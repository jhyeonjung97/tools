for model in gpr gbr rf lr xgb lgb; do
    for numb in 0.7 0.8 0.9 1.0; do
    
python ~/bin/tools/tetra/bulk_pred.py --model lr --corr_threshold 0.7
python ~/bin/tools/tetra/bulk_pred.py --model lr --corr_threshold 0.8
python ~/bin/tools/tetra/bulk_pred.py --model lr --corr_threshold 0.9
python ~/bin/tools/tetra/bulk_pred.py --model lr

python ~/bin/tools/tetra/bulk_pred.py --model rf --corr_threshold 0.7
python ~/bin/tools/tetra/bulk_pred.py --model rf --corr_threshold 0.8
python ~/bin/tools/tetra/bulk_pred.py --model rf --corr_threshold 0.9
python ~/bin/tools/tetra/bulk_pred.py --model rf

python ~/bin/tools/tetra/bulk_pred.py --model gpr --corr_threshold 0.7
python ~/bin/tools/tetra/bulk_pred.py --model gpr --corr_threshold 0.8
python ~/bin/tools/tetra/bulk_pred.py --model gpr --corr_threshold 0.9
python ~/bin/tools/tetra/bulk_pred.py --model gpr

python ~/bin/tools/tetra/bulk_pred.py --model gbr --corr_threshold 0.7
python ~/bin/tools/tetra/bulk_pred.py --model gbr --corr_threshold 0.8
python ~/bin/tools/tetra/bulk_pred.py --model gbr --corr_threshold 0.9
python ~/bin/tools/tetra/bulk_pred.py --model gbr

python ~/bin/tools/tetra/bulk_pred.py --model xgb --corr_threshold 0.7
python ~/bin/tools/tetra/bulk_pred.py --model xgb --corr_threshold 0.8
python ~/bin/tools/tetra/bulk_pred.py --model xgb --corr_threshold 0.9
python ~/bin/tools/tetra/bulk_pred.py --model xgb

python ~/bin/tools/tetra/bulk_pred.py --model lgb --corr_threshold 0.7
python ~/bin/tools/tetra/bulk_pred.py --model lgb --corr_threshold 0.8
python ~/bin/tools/tetra/bulk_pred.py --model lgb --corr_threshold 0.9
python ~/bin/tools/tetra/bulk_pred.py --model lgb