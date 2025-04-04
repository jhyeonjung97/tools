#!/bin/bash

# 현재 날짜와 시간으로 로그 파일명 생성
LOG_FILE="coord_pred_$(date +%Y%m%d_%H%M%S).log"

# 스크립트 실행 및 출력을 화면과 로그 파일에 동시에 저장
{
    echo "Starting coordination prediction analysis at $(date)"
    echo "----------------------------------------"
    
    for model in gb lr svm rf; do
        echo "Running model: $model"
        python ~/bin/tools/tetra/coord_pred.py --model $model 2>&1
        # for row in 3d 4d 5d '4d 5d'; do
        #     echo "Running model: $model, row: $row"
        #     python ~/bin/tools/tetra/coord_pred.py --model $model --row $row 2>&1
        # done
    done
    
    echo "----------------------------------------"
    echo "Finished coordination prediction analysis at $(date)"
} 2>&1 | tee $LOG_FILE
