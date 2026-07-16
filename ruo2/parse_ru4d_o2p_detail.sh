#!/bin/bash
# Break down ICOHP(Ru4d-O2p) per Ru-O bond and per Ru atom (spin1+spin2 total),
# in the current directory's ICOHPLIST.lobster.
# Some ICOHPLIST.lobster files list spin1 and spin2 as two columns on the same
# line (sum both columns); others list a full spin1 block followed by a full
# spin2 block, each with its own "for spin N" header (sum every matching line
# across both blocks, since each line already carries one spin's value).
# Atom labels are shifted -1 (Ru1->Ru0, O9->O8, ...) to match ASE's
# 0-indexed atom numbering, since LOBSTER numbers atoms starting at 1.
if [ ! -e ICOHPLIST.lobster ]; then
    exit 0
fi

if grep -q "for spin  2" ICOHPLIST.lobster; then
    two_col=0
else
    two_col=1
fi

awk -v OFS='\t' -v two_col="$two_col" '
function ase_label(label,   n) {
    n = label
    gsub(/[A-Za-z]/, "", n)
    gsub(n"$", "", label)
    return label (n - 1)
}
function bond_icohp(   v) {
    v = 0
    if (two_col) {
        if ($(NF-1) ~ /^-?[0-9.]+$/) v += $(NF-1)
        if ($NF ~ /^-?[0-9.]+$/) v += $NF
    } else {
        if ($NF ~ /^-?[0-9.]+$/) v += $NF
    }
    return v
}
($2 ~ /^Ru[0-9]+_4d/ && $3 ~ /^O[0-9]+_2p/) {
    ru = $2; sub(/_4d.*/, "", ru); ru = ase_label(ru)
    o  = $3; sub(/_2p.*/, "", o); o = ase_label(o)
    key = $1 SUBSEP ru SUBSEP o
    if (!(key in dist)) dist[key] = $4
    sum[key] += bond_icohp()
    next
}
($3 ~ /^Ru[0-9]+_4d/ && $2 ~ /^O[0-9]+_2p/) {
    ru = $3; sub(/_4d.*/, "", ru); ru = ase_label(ru)
    o  = $2; sub(/_2p.*/, "", o); o = ase_label(o)
    key = $1 SUBSEP ru SUBSEP o
    if (!(key in dist)) dist[key] = $4
    sum[key] += bond_icohp()
    next
}
END {
    print "COHP_num","Ru_atom","O_atom","distance","ICOHP_sum" > "icohp_Ru4d_O2p_perbond.tsv"
    for (k in sum) {
        split(k, a, SUBSEP)
        print a[1], a[2], a[3], dist[k], sprintf("%.2f", sum[k]) >> "icohp_Ru4d_O2p_perbond.tsv"
        peratom[a[2]] += sum[k]
        nbond[a[2]]++
    }
    print "Ru_atom","ICOHP_sum","n_bonds" > "icohp_Ru4d_O2p_peratom.tsv"
    for (r in peratom) {
        print r, sprintf("%.2f", peratom[r]), nbond[r] >> "icohp_Ru4d_O2p_peratom.tsv"
    }
}
' ICOHPLIST.lobster

{ head -1 icohp_Ru4d_O2p_perbond.tsv; tail -n +2 icohp_Ru4d_O2p_perbond.tsv | sort -t$'\t' -k2,2V -k3,3V; } > icohp_Ru4d_O2p_perbond.tsv.tmp && mv icohp_Ru4d_O2p_perbond.tsv.tmp icohp_Ru4d_O2p_perbond.tsv
{ head -1 icohp_Ru4d_O2p_peratom.tsv; tail -n +2 icohp_Ru4d_O2p_peratom.tsv | sort -t$'\t' -k1,1V; } > icohp_Ru4d_O2p_peratom.tsv.tmp && mv icohp_Ru4d_O2p_peratom.tsv.tmp icohp_Ru4d_O2p_peratom.tsv
