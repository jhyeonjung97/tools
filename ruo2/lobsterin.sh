for dir in */; do
    cd $dir
    dir=${dir%/}
    IFS='_' read -r A B <<< "$dir"
    sed -i -e "s/X/$B/g" lobsterin
    cd $dir_now
done