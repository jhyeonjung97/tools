for model in lr gpr gbr rf; do
    for Y in form coh; do
    python ~/bin/tools/tetra/bulk_pred.py --model $model --Y $Y
        for row in 3d 4d 5d '4d 5d'; do
            python ~/bin/tools/tetra/bulk_pred.py --model $model --Y $Y --row $row
        done
    done
done