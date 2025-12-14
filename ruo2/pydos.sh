dir1="/Users/jiuy97/Desktop/3_RuO2/1_RuO2_Ueff/0_bulk/between_2_and_3/5_"
dir2="/Users/jiuy97/Desktop/3_RuO2/3_ReRuO2_OER/0_bulk"

# Save current directory
original_dir=$(pwd)

cd $dir1
fermi1=$(grep "E-fermi" OUTCAR 2>/dev/null | tail -n 1 | awk '{print $3}')

cd $dir2
fermi2=$(grep "E-fermi" OUTCAR 2>/dev/null | tail -n 1 | awk '{print $3}')

# Use fermi1 as default if fermi2 is empty
if [ -z "$fermi2" ]; then
    fermi2=$fermi1
fi

# Output file path (save in original directory)
output_file="$original_dir/dos-bulk.pdf"

python ~/bin/tools/ruo2/pydos-bulk.py \
-d "$dir2" --dir2 "$dir1" \
-z "$fermi2" --zero2 "$fermi1" \
-y -12 12 --notot \
--elem2='Ru' --spd2='d' --lc2 'C2' \
--elem2='O' --spd2='p' --lc2 'C1' \
--elem='Re' --spd='d' --lc 'C0' \
--elem='Ru' --spd='d' --lc 'C2' \
--elem='O' --spd='p' --lc 'C1' \
-o "$output_file"

echo "Output saved to: $output_file"