awk '$2 !~ /_/ && $3 !~ /_/' ICOHPLIST.lobster > icohp.txt

awk '$2 ~ /2s/ && $3 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_4s.txt
awk '$2 ~ /2s/ && $3 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_5s.txt
awk '$2 ~ /2s/ && $3 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_6s.txt

awk '$2 ~ /2p/ && $3 ~ /4s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_4s.txt
awk '$2 ~ /2p/ && $3 ~ /5s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_5s.txt
awk '$2 ~ /2p/ && $3 ~ /6s/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_6s.txt

awk '$2 ~ /2s/ && $3 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_3d.txt
awk '$2 ~ /2s/ && $3 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_4d.txt
awk '$2 ~ /2s/ && $3 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2s_5d.txt

awk '$2 ~ /2p/ && $3 ~ /3d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_3d.txt
awk '$2 ~ /2p/ && $3 ~ /4d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_4d.txt
awk '$2 ~ /2p/ && $3 ~ /5d/ && $NF ~ /^-?[0-9.]+$/ {sum += $NF} END {print sum}' ICOHPLIST.lobster > icohp_2p_5d.txt

for file in icohp*.txt; do
    [ -e "$file" ] || continue
    if [ ! -s "$file" ]; then
        rm "$file"
    fi
done