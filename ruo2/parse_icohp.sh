# awk '$2 !~ /_/ && $3 !~ /_/' ICOHPLIST.lobster > icohp.txt

if [ -e icohp_sum.txt ]; then
    rm icohp_sum.txt
fi
if [ -e icohp_dopant.txt ]; then
    rm icohp_dopant.txt
fi

awk '$2 ~ /2p/ && $3 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$2 ~ /2p/ && $3 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$2 ~ /2p/ && $3 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$3 ~ /2p/ && $2 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$3 ~ /2p/ && $2 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$3 ~ /2p/ && $2 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$2 ~ /2p/ && $3 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$2 ~ /2p/ && $3 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$2 ~ /2p/ && $3 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$3 ~ /2p/ && $2 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$3 ~ /2p/ && $2 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt
awk '$3 ~ /2p/ && $2 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_sum.txt

awk '$2 ~ /2p/ && $3 ~ /4s/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$2 ~ /2p/ && $3 ~ /5s/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$2 ~ /2p/ && $3 ~ /6s/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$3 ~ /2p/ && $2 ~ /4s/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$3 ~ /2p/ && $2 ~ /5s/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$3 ~ /2p/ && $2 ~ /6s/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$2 ~ /2p/ && $3 ~ /3d/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$2 ~ /2p/ && $3 ~ /4d/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$2 ~ /2p/ && $3 ~ /5d/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$3 ~ /2p/ && $2 ~ /3d/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$3 ~ /2p/ && $2 ~ /4d/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt
awk '$3 ~ /2p/ && $2 ~ /5d/ && $2 !~ /Ru/ && $3 !~ /Ru/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_dopant.txt

# awk '$2 ~ /2s/ && $3 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_4s.txt
# awk '$2 ~ /2s/ && $3 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_5s.txt
# awk '$2 ~ /2s/ && $3 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_6s.txt

# awk '$2 ~ /2p/ && $3 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_4s.txt
# awk '$2 ~ /2p/ && $3 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_5s.txt
# awk '$2 ~ /2p/ && $3 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_6s.txt

# awk '$2 ~ /2s/ && $3 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_3d.txt
# awk '$2 ~ /2s/ && $3 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_4d.txt
# awk '$2 ~ /2s/ && $3 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_5d.txt

# awk '$2 ~ /2p/ && $3 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_3d.txt
# awk '$2 ~ /2p/ && $3 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_4d.txt
# awk '$2 ~ /2p/ && $3 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_5d.txt

# awk '$3 ~ /2s/ && $2 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2s_4s.txt
# awk '$3 ~ /2s/ && $2 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2s_5s.txt
# awk '$3 ~ /2s/ && $2 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2s_6s.txt

# awk '$3 ~ /2p/ && $2 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2p_4s.txt
# awk '$3 ~ /2p/ && $2 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2p_5s.txt
# awk '$3 ~ /2p/ && $2 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2p_6s.txt

# awk '$3 ~ /2s/ && $2 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2s_3d.txt
# awk '$3 ~ /2s/ && $2 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2s_4d.txt
# awk '$3 ~ /2s/ && $2 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2s_5d.txt

# awk '$3 ~ /2p/ && $2 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2p_3d.txt
# awk '$3 ~ /2p/ && $2 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2p_4d.txt
# awk '$3 ~ /2p/ && $2 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster >> icohp_2p_5d.txt

for file in icohp*.txt; do
    [ -e "$file" ] || continue
    sed -i '/^\s*$/d' "$file"
done

for file in icohp*.txt; do
    [ -e "$file" ] || continue
    if [ "$(stat -c%s "$file")" -le 2 ]; then
        rm "$file"
    fi
done