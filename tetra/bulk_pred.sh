#!/bin/bash

for model in gpr gbr rf lr; do
    for numb in 0 7 8 9; do
        python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold7_old0 --X chg dipole density ion1 volume ICOOPn madelung mag Hfus n_bond grosspop
        python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold8_old0 --X chg pauling ion2 dipole density numb ion1 Vatom Natom ICOOPn madelung l_bond mag Hfus ion3 n_bond grosspop
        python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold9_old0 --X chg pauling ion12 dipole density numb Rcoval ion1 Vatom Natom Rvdw volume ICOOPn ICOBIm madelung l_bond mag Hfus ion3 Hevap n_bond grosspop
        python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold0_old0
    done
done

for model in lr; do
    python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold7_old0 --X chg dipole density ion1 volume ICOOPn madelung mag Hfus n_bond grosspop
    python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold8_old0 --X chg pauling ion2 dipole density numb ion1 Vatom Natom ICOOPn madelung l_bond mag Hfus ion3 n_bond grosspop
    python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold9_old0 --X chg pauling ion12 dipole density numb Rcoval ion1 Vatom Natom Rvdw volume ICOOPn ICOBIm madelung l_bond mag Hfus ion3 Hevap n_bond grosspop
    python ~/bin/tools/tetra/bulk_pred0.py --model $model --Y form --output threshold0_old0
done

python ~/bin/tools/tetra/bulk_pred0.py --model gpr --output threshold0_old0
python ~/bin/tools/tetra/bulk_pred7.py --model gpr --output threshold0_old7
python ~/bin/tools/tetra/bulk_pred8.py --model gpr --output threshold0_old8
python ~/bin/tools/tetra/bulk_pred9.py --model gpr --output threshold0_old9

# # 현재 날짜와 시간으로 로그 파일명 생성
# LOG_FILE="bulk_pred_$(date +%Y%m%d_%H%M%S).log"

# # 스크립트 실행 및 출력을 화면과 로그 파일에 동시에 저장
# {
#     echo "Starting bulk prediction analysis at $(date)"
#     echo "----------------------------------------"
    
#     for model in gpr gbr rf lr; do
#         for Y in form coh; do
#             echo "Running model: $model target: $Y"
#             python ~/bin/tools/tetra/bulk_pred.py --model $model --Y $Y 2>&1
#             for row in 3d 4d 5d '4d 5d'; do
#                 echo "Running model: $model target: $Y row: $row"
#                 python ~/bin/tools/tetra/bulk_pred.py --model $model --Y $Y --row $row 2>&1
#             done
#         done
#     done
    
#     echo "----------------------------------------"
#     echo "Finished bulk prediction analysis at $(date)"
# } 2>&1 | tee $LOG_FILE