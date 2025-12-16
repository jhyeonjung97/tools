#!/bin/bash

for file in */; do
    cd "$file"
    
    # 현재 경로를 '/'로 분할하여 배열로 저장
    IFS='/' read -r -a path_components <<< "$PWD"
    
    # 배열의 길이를 확인하고 적절한 인덱스 사용
    array_length=${#path_components[@]}
    if [ $array_length -ge 2 ]; then
        dir=${path_components[$((array_length-1))]}
        # numb=${path_components[$((array_length-1))]}
        
        name=${dir#*_}
        suffix=${numb%_}
        
        echo "Directory: $dir"
        echo "Number: $numb"
        echo "Name: $name"
        echo "Suffix: $suffix"
        echo "---"
    fi

    # 원래 디렉토리로 돌아가기
    cd - > /dev/null
    cp ./${file}restart.json ./${name}.json

done

rm Ru_top.json
