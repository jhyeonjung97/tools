for dir in */; do
    IFS='_' read -r A B <<< "$dir"
    echo $B
    # sed -i -e "s/X/$B/g" $dir/lobsterin
done