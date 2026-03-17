dirs=('1_Hf' '2_Ta' '3_W' '4_Re' '5_Os')
xs=(a2 a4 a6 a8)
ys=(a b c d)
for i in {0..4}; do
    for j in {0..3}; do
        k=$((i+1))
        cp ${dirs[$i]}/${xs[$j]}.json ${ys[$j]}${k}.json
    done
done