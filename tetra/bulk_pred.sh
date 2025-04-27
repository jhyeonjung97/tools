#!/bin/bash

for model in rf lr; do
    for numb in 0.5 0.6 0.7 0.8 0.9 1.0; do
        python ~/bin/tools/tetra/bulk_pred0.py --model $model --corr_threshold $numb
        python ~/bin/tools/tetra/bulk_pred_cfse.py --model $model --corr_threshold $numb
    done
done

for model in gpr gbr xgb lgb; do
    for numb in $(seq 51 99 | awk '{printf "%.2f\n", $1/100}'); do
        python ~/bin/tools/tetra/bulk_pred0.py --model $model --corr_threshold $numb
        python ~/bin/tools/tetra/bulk_pred_cfse.py --model $model --corr_threshold $numb
    done
done

for numb in 0.5 0.6 0.7 0.8 0.9 1.0; do
    python ~/bin/tools/tetra/bulk_pred0.py --model lr --corr_threshold $numb
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model lr --corr_threshold $numb
done

for numb in 0.5 0.6 0.7 0.8 0.9 1.0; do
    python ~/bin/tools/tetra/bulk_pred0.py --model rf --corr_threshold $numb
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model rf --corr_threshold $numb
done

for numb in $(seq 31 39 | awk '{printf "%.2f\n", $1/100}'); do
    python ~/bin/tools/tetra/bulk_pred0.py --model gpr --corr_threshold $numb
done

for numb in $(seq 41 49 | awk '{printf "%.2f\n", $1/100}'); do
    python ~/bin/tools/tetra/bulk_pred0.py --model gpr --corr_threshold $numb
done

for numb in $(seq 31 39 | awk '{printf "%.2f\n", $1/100}'); do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --corr_threshold $numb
done

for numb in $(seq 41 49 | awk '{printf "%.2f\n", $1/100}'); do
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gpr --corr_threshold $numb
done

for numb in $(seq 31 49 | awk '{printf "%.2f\n", $1/100}'); do
    python ~/bin/tools/tetra/bulk_pred0.py --model gbr --corr_threshold $numb
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model gbr --corr_threshold $numb
done

for numb in $(seq 31 49 | awk '{printf "%.2f\n", $1/100}'); do
    python ~/bin/tools/tetra/bulk_pred0.py --model xgb --corr_threshold $numb
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model xgb --corr_threshold $numb
done

for numb in $(seq 31 49 | awk '{printf "%.2f\n", $1/100}'); do
    python ~/bin/tools/tetra/bulk_pred0.py --model lgb --corr_threshold $numb
    python ~/bin/tools/tetra/bulk_pred_cfse.py --model lgb --corr_threshold $numb
done