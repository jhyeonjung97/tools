python ~/bin/verve/lr.py -o Ef \
--Y concat_norm_formation.tsv \
-i \
ICOBI_per_metal \
ICOBI_per_bond \
chg \
en_pauling \
density \
vdw_radius \
fusion_heat \
--X \
concat_ICOBI.tsv \
concat_norm_ICOBI.tsv \
concat_chg.tsv \
concat_en_pauling.tsv \
concat_density.tsv \
concat_vdw_radius.tsv \
concat_fusion_heat.tsv