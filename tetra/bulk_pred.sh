#!/bin/bash

# 현재 날짜와 시간으로 로그 파일명 생성
LOG_FILE="bulk_pred_$(date +%Y%m%d_%H%M%S).log"

# 스크립트 실행 및 출력을 화면과 로그 파일에 동시에 저장
{
    echo "Starting bulk prediction analysis at $(date)"
    echo "----------------------------------------"
    
    for model in gpr gbr rf lr; do
        for Y in form coh; do
            echo "Running model: $model, target: $Y"
            python ~/bin/tools/tetra/bulk_pred.py --model $model --Y $Y 2>&1
            for row in 3d 4d 5d '4d 5d'; do
                echo "Running model: $model, target: $Y, row: $row"
                python ~/bin/tools/tetra/bulk_pred.py --model $model --Y $Y --row $row 2>&1
            done
        done
    done
    
    echo "----------------------------------------"
    echo "Finished bulk prediction analysis at $(date)"
} 2>&1 | tee $LOG_FILE