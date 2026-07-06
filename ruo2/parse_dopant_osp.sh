#!/bin/bash
# Sum ICOHP (spin1+spin2 total) for dopant(shell)-O(2s+2p) bonds in the
# current directory's ICOHPLIST.lobster.
# Usage: parse_dopant_osp.sh <ElementSymbol> <shell e.g. 3d|4d|5d>
# Some ICOHPLIST.lobster files list spin1 and spin2 as two columns on the same
# line (sum both columns); others list a full spin1 block followed by a full
# spin2 block, each with its own "for spin N" header (sum every matching line
# across both blocks, since each line already carries one spin's value).
el="$1"
shell="$2"
if [ -z "$el" ] || [ -z "$shell" ] || [ ! -e ICOHPLIST.lobster ]; then
    echo ""
    exit 0
fi

if grep -q "for spin  2" ICOHPLIST.lobster; then
    two_col=0
else
    two_col=1
fi

awk -v el="$el" -v shell="$shell" -v two_col="$two_col" '
(($2 ~ "^"el"[0-9]+_"shell) && ($3 ~ /^O[0-9]+_2[sp]/)) ||
(($3 ~ "^"el"[0-9]+_"shell) && ($2 ~ /^O[0-9]+_2[sp]/)) {
    if (two_col) {
        if ($(NF-1) ~ /^-?[0-9.]+$/) sum += $(NF-1)
        if ($NF ~ /^-?[0-9.]+$/) sum += $NF
    } else {
        if ($NF ~ /^-?[0-9.]+$/) sum += $NF
    }
}
END { printf "%.2f\n", sum+0 }
' ICOHPLIST.lobster
