#!/bin/bash

for model in gpr gbr rf lr xgb lgb; do
    for numb in 0.7 0.8 0.9 1.0; do
        python ~/bin/tools/tetra/bulk_pred.py --model $model --corr_threshold $numb
    done
done

/opt/homebrew/bin/rsync -avz '/Users/jiuy97/Desktop/7_V_bulk/figures' '/Users/jiuy97/Library/CloudStorage/GoogleDrive-jiuy97@stanford.edu/My Drive/Tetrahedral_oxides_ML/Figures'