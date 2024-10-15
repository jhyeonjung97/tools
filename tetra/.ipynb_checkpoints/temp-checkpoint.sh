### section5-2

cd /pscratch/sd/j/jiuy97/3_V_bulk/figures
python ~/bin/tools/tetra/concat.py -o energy --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_energy.tsv
python ~/bin/tools/tetra/concat.py -o norm_energy --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_norm_energy.tsv
python ~/bin/tools/tetra/concat.py -o norm_formation --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_norm_formation.tsv
python ~/bin/tools/tetra/concat.py -o norm_cohesive --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_norm_cohesive.tsv
python ~/bin/tools/tetra/concat.py -o ICOHP --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_ICOHP.tsv
python ~/bin/tools/tetra/concat.py -o norm_ICOHP --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_norm_ICOHP.tsv
python ~/bin/tools/tetra/concat.py -o ICOBI --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_ICOBI.tsv
python ~/bin/tools/tetra/concat.py -o norm_ICOBI --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_norm_ICOBI.tsv
python ~/bin/tools/tetra/concat.py -o bond --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_bond.tsv
python ~/bin/tools/tetra/concat.py -o norm_volume --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_norm_volume.tsv
python ~/bin/tools/tetra/concat.py -o chg --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_chg.tsv
python ~/bin/tools/tetra/concat.py -o mag --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_mag_M.tsv
python ~/bin/tools/tetra/concat.py -o norm_MadelungL --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_norm_Madelung_L.tsv
python ~/bin/tools/tetra/concat.py -o GrossPopulationL --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_GP_L_M.tsv
python ~/bin/tools/tetra/concat.py -o redoxP --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_redoxP.tsv
python ~/bin/tools/tetra/concat.py -o redoxP_clean --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_redoxP_clean.tsv
python ~/bin/tools/tetra/rel2octa.py concat_bond.tsv
python ~/bin/tools/tetra/rel2octa.py concat_chg.tsv
python ~/bin/tools/tetra/rel2octa.py concat_energy.tsv
python ~/bin/tools/tetra/rel2octa.py concat_GrossPopulationL.tsv
python ~/bin/tools/tetra/rel2octa.py concat_ICOBI.tsv
python ~/bin/tools/tetra/rel2octa.py concat_ICOHP.tsv
python ~/bin/tools/tetra/rel2octa.py concat_mag.tsv
python ~/bin/tools/tetra/rel2octa.py concat_norm_energy.tsv
python ~/bin/tools/tetra/rel2octa.py concat_norm_formation.tsv
python ~/bin/tools/tetra/rel2octa.py concat_norm_cohesive.tsv
python ~/bin/tools/tetra/rel2octa.py concat_norm_ICOHP.tsv
python ~/bin/tools/tetra/rel2octa.py concat_norm_ICOBI.tsv
python ~/bin/tools/tetra/rel2octa.py concat_norm_MadelungL.tsv
python ~/bin/tools/tetra/rel2octa.py concat_norm_volume.tsv
mv concat*rel.tsv rel5/
mv concat*.tsv rel6/

### section5-3

cd /pscratch/sd/j/jiuy97/3_V_bulk/figures/rel5
python ~/bin/tools/tetra/mendeleev2tsv.py -n 5 -p \
group_id atomic_number atomic_volume  \
boiling_point melting_point \
mass density dipole_polarizability en_pauling \
covalent_radius metallic_radius vdw_radius \
evaporation_heat fusion_heat  heat_of_formation \
ionenergies[1] ionenergies[2] ionenergies[3]
python ~/bin/tools/tetra/operator.py -o + \
-x concat_ionenergies_1.tsv \
-y concat_ionenergies_2.tsv \
-z concat_ionenergies_12.tsv
python ~/bin/tools/tetra/operator.py -o + \
-x concat_evaporation_heat.tsv \
-y concat_fusion_heat.tsv \
-z concat_sublimation_heat.tsv
declare -A files_A
files_A[coord]="merged_coord.tsv"
files_A[redoxP]="merged_redoxP.tsv"
files_A[redoxP_clean]="merged_redoxP_clean.tsv"
for key in "${!files_A[@]}"; do
    python ~/bin/tools/tetra/concat.py -o $key --X \
    /pscratch/sd/j/jiuy97/3_V_bulk/1_Tetrahedral_WZ/${files_A[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/2_Tetrahedral_ZB/${files_A[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/3_SquarePlanar_TN/${files_A[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/4_SquarePlanar_PD/${files_A[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/5_SquarePlanar_NB/${files_A[$key]} \
done
declare -A files_B
files_B[element]="merged_element.tsv"
files_B[row]="merged_row.tsv"
for key in "${!files_B[@]}"; do
    python ~/bin/tools/tetra/concat.py -o $key --X \
    /pscratch/sd/j/jiuy97/3_V_bulk/metal/${files_B[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/metal/${files_B[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/metal/${files_B[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/metal/${files_B[$key]} \
    /pscratch/sd/j/jiuy97/3_V_bulk/metal/${files_B[$key]}
done
for file in ~/bin/tools/tetra/png_rel/lr*.sh; do
    sh $file
done

### section5-4

cd /pscratch/sd/j/jiuy97/3_V_bulk/figures/rel6
python ~/bin/tools/tetra/mendeleev2tsv.py -n 6 -p \
group_id atomic_number atomic_volume  \
boiling_point melting_point \
mass density dipole_polarizability en_pauling \
covalent_radius metallic_radius vdw_radius \
evaporation_heat fusion_heat  heat_of_formation \
ionenergies[1] ionenergies[2] ionenergies[3]
python ~/bin/tools/tetra/operator.py -o + \
-x concat_ionenergies_1.tsv \
-y concat_ionenergies_2.tsv \
-z concat_ionenergies_12.tsv
python ~/bin/tools/tetra/operator.py -o + \
-x concat_evaporation_heat.tsv \
-y concat_fusion_heat.tsv \
-z concat_sublimation_heat.tsv
python ~/bin/tools/tetra/concat.py -o coord --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_coord.tsv
python ~/bin/tools/tetra/concat.py -o element --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_element.tsv
python ~/bin/tools/tetra/concat.py -o row --X /pscratch/sd/j/jiuy97/3_V_bulk/*_*_*/merged_row.tsv
for file in ~/bin/tools/tetra/png/lr*.sh; do
    sh $file
done