for dir in /Users/jiuy97/Desktop/8_V_slab/*_*_*/*/*_*
do
    if [ -f "$dir/unmatched" ]; then
        rm -r $dir/*.*
    fi
done

for dir in /Users/jiuy97/Desktop/8_V_slab/*_*_*/*/*_*
do
    if [ ! -f "$dir/DONE" ] && [ -f "$dir/final_with_calculator.json" ]; then
        echo "$dir"
    fi
done
