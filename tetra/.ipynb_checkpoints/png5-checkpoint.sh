#!/bin/bash

### section1

# for dir in /Users/hailey/Desktop/3_V_bulk/metal/*/
# do
#     cd $dir
#     python ~/bin/tools/tetra/energy.py --save -x "Metal (M)" -y "Total energy (eV/M)" -n m
# done
# cd /Users/hailey/Desktop/3_V_bulk/metal
# python ~/bin/tools/tetra/tsv.py -l  3d 4d 5d -x "Metal (MO)" -y "Total energy (eV)" -o norm_energy *d/energy_norm_energy.tsv
# cd /Users/hailey/Desktop/3_V_bulk/oxide/0_min
# python ~/bin/tools/tetra/energy.py --save -x "Metal (MxOy)" -y "Total energy (eV/M)" -n m
# cd /Users/hailey/Desktop/3_V_bulk/6_Octahedral_RS
# python ~/bin/tools/tetra/mendeleev2tsv.py -p \
# group_id atomic_number atomic_volume  \
# boiling_point melting_point \
# mass density dipole_polarizability en_pauling \
# covalent_radius metallic_radius vdw_radius \
# evaporation_heat fusion_heat heat_of_formation \
# ionenergies[1] ionenergies[2] ionenergies[3]
# python ~/bin/tools/tetra/operator.py -o + \
# -x mendeleev_ionenergies_1.tsv \
# -y mendeleev_ionenergies_2.tsv \
# -z mendeleev_ionenergies_12.tsv
# python ~/bin/tools/tetra/operator.py -o + \
# -x mendeleev_evaporation_heat.tsv \
# -y mendeleev_fusion_heat.tsv \
# -z mendeleev_sublimation_heat.tsv
# if [[ ${here} == 'nersc' ]]; then
#     sed -i 's/\x0//g' *.tsv
# fi

### section2

for dir in /Users/hailey/Desktop/3_V_bulk/6_Octahedral_RS/*/
do
    cd $dir
    if [[ $dir == *'Tetrahedral'* ]]; then
        n=4; #python ~/bin/tools/tetra/energy.py --save -p hexa -x "Metal (MO)" -y "Hexagonal ratio [c/a]"
    elif [[ $dir == *'Pyramidal'* ]] || [[ $dir == *'Tetrahedral'* ]] || [[ $dir == *'SquarePlanar'* ]]; then
        n=4; #python ~/bin/tools/tetra/energy.py --save -p hexa -x "Metal (MO)" -y "Square prism ratio [c/a]"
    elif [[ $dir == *'Octahedral'* ]]; then
        n=6
    fi
    python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV)"
    python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV/MO)" -n m
    python ~/bin/tools/tetra/energy.py --save -p bond  -x "Metal (MO)" -y "Bond length (A)"
    python ~/bin/tools/tetra/energy.py --save -p bond -x "Metal (MO)" -y "Bond length (A/M-O)" -n $n
    python ~/bin/tools/tetra/energy.py --save -p volume -x "Metal (MO)" -y "Volume (A^3/MO)" -n m
    python ~/bin/tools/tetra/energy.py --save -p chg -e M  -x "Metal (MO)" -y "Bader charge (e-)"
    python ~/bin/tools/tetra/energy.py --save -p mag -e M -x "Metal (MO)" -y "|Magnetization|"
    python ~/bin/tools/tetra/energy.py --save -p ICOHP -x "Metal (MO)" -y "ICOHP (eV/MO)"
    python ~/bin/tools/tetra/energy.py --save -p ICOHP -x "Metal (MO)" -y "ICOHP (eV/M-O)" -n $n
    python ~/bin/tools/tetra/energy.py --save -p ICOBI -x "Metal (MO)" -y "ICOBI (/MO)"
    python ~/bin/tools/tetra/energy.py --save -p ICOBI -x "Metal (MO)" -y "ICOBI (eV/M-O)" -n $n
    python ~/bin/tools/tetra/energy.py --save -p GP_L -e M  -x "Metal (MO)" -y "Gross population (Loewdin)"
    python ~/bin/tools/tetra/energy.py --save -p Madelung_L -x "Metal (MO)" -y "Madelugn energy (Loewdin, eV/MO)" -n m
    python ~/bin/tools/tetra/formation_energy.py
    python ~/bin/tools/tetra/cohesive_energy.py
    if [[ ${here} == 'nersc' ]]; then
        sed -i 's/\x0//g' *.tsv
    fi
done

### section3

for dir in /Users/hailey/Desktop/3_V_bulk/6_Octahedral_RS/
do
    cd $dir
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Total energy (eV)" \
    -o energy 3d/energy_energy.tsv 4d/energy_energy.tsv 5d/energy_energy.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Total energy (eV/MO)" \
    -o norm_energy 3d/energy_norm_energy.tsv 4d/energy_norm_energy.tsv 5d/energy_norm_energy.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Bond length (A/M-O)" \
    -o bond 3d/energy_bond.tsv 4d/energy_bond.tsv 5d/energy_bond.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Volume (A^3/MO)" \
    -o norm_volume 3d/energy_norm_volume.tsv 4d/energy_norm_volume.tsv 5d/energy_norm_volume.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Bader charge (e-)" \
    -o chg 3d/energy_chg_M.tsv 4d/energy_chg_M.tsv 5d/energy_chg_M.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "|Magnetization|" \
    -o mag_M 3d/energy_mag_M.tsv 4d/energy_mag_M.tsv 5d/energy_mag_M.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOHP (eV/MO)" \
    -o ICOHP 3d/energy_ICOHP.tsv 4d/energy_ICOHP.tsv 5d/energy_ICOHP.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOHP (eV/M-O)" \
    -o norm_ICOHP 3d/energy_norm_ICOHP.tsv 4d/energy_norm_ICOHP.tsv 5d/energy_norm_ICOHP.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOBI (/MO)" \
    -o ICOBI 3d/energy_ICOBI.tsv 4d/energy_ICOBI.tsv 5d/energy_ICOBI.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOBI (/M-O)" \
    -o norm_ICOBI 3d/energy_norm_ICOBI.tsv 4d/energy_norm_ICOBI.tsv 5d/energy_norm_ICOBI.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Gross population (Loewdin)" \
    -o GP_L_M 3d/energy_GP_Loewdin_M.tsv 4d/energy_GP_Loewdin_M.tsv 5d/energy_GP_Loewdin_M.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Madelung energy (Loewdin, eV/MO)" \
    -o norm_Madelung_L 3d/energy_norm_Madelung_Loewdin.tsv 4d/energy_norm_Madelung_Loewdin.tsv 5d/energy_norm_Madelung_Loewdin.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Formation energy (eV/MO)" \
    -o norm_formation 3d/energy_norm_formation.tsv 4d/energy_norm_formation.tsv 5d/energy_norm_formation.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Cohesive energy (eV/MO)" \
    -o norm_cohesive 3d/energy_norm_cohesive.tsv 4d/energy_norm_cohesive.tsv 5d/energy_norm_cohesive.tsv
    python ~/bin/tools/tetra/tsv.py -l 3d fm -x "Metal (MO)" -y "Formation energy (eV/MO)" \
    -o AFMvsFM 3d/energy_norm_formation.tsv fm/energy_norm_formation.tsv
done

### section4-1

# cd /Users/hailey/Desktop/3_V_bulk/6_Octahedral_RS
# python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Standard reduction potential (V)" -o redoxP \
# ../oxide/energy_redoxP_3d.tsv ../oxide/energy_redoxP_4d.tsv ../oxide/energy_redoxP_5d.tsv
# python ~/bin/tools/tetra/tsv.py -l 3d 4d 5d -x "Metal (MO)" -y "Standard reduction potential (V)" -o redoxP_clean \
# ../oxide/energy_redoxP_clean_3d.tsv ../oxide/energy_redoxP_clean_4d.tsv ../oxide/energy_redoxP_clean_5d.tsv

### section4-2

# cd /Users/hailey/Desktop/3_V_bulk
# sh ~/bin/verve/spread.sh 6_Octahedral_RS/merged_redoxP.tsv
# sh ~/bin/verve/spread.sh 6_Octahedral_RS/merged_redoxP_clean.tsv
# sh ~/bin/verve/spread.sh oxide/merged_coord.tsv
# sh ~/bin/verve/spread.sh metal/merged_element.tsv
# sh ~/bin/verve/spread.sh metal/merged_row.tsv
# sed -i -e 's/CN/WZ/g' 1_Tetrahedral_WZ/merged_coord.tsv
# sed -i -e 's/CN/ZB/g' 2_Tetrahedral_ZB/merged_coord.tsv
# sed -i -e 's/CN/TN/g' 3_SquarePlanar_TN/merged_coord.tsv
# sed -i -e 's/CN/PD/g' 4_SquarePlanar_PD/merged_coord.tsv
# sed -i -e 's/CN/NB/g' 5_SquarePlanar_NB/merged_coord.tsv
# sed -i -e 's/CN/RS/g' 6_Octahedral_RS/merged_coord.tsv
# sed -i -e 's/CN/LT/g' 7_Pyramidal_LT/merged_coord.tsv
# sed -i -e 's/CN/AQ/g' 8_Tetrahedral_AQ/merged_coord.tsv
# sed -i -e 's/CN/AU/g' 9_SquarePlanar_AU/merged_coord.tsv

### section5-1

cd /Users/hailey/Desktop/3_V_bulk/figures
for row in 3d 4d 5d
do
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Total energy (eV)" \
    -o energy_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_energy.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Total energy (eV/MO)" \
    -o norm_energy_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_norm_energy.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Formation energy (eV/MO)" \
    -o norm_formation_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_norm_formation.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Cohesive energy (eV/MO)" \
    -o norm_cohesive_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_norm_cohesive.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Bond length (A/M-O)" \
    -o bond_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_bond.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Volume (A^3/MO)" \
    -o norm_volume_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_norm_volume.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Bader charge (e-)" \
    -o chg_M_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_chg_M.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "|Magnetization|" \
    -o mag_M_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_mag_M.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "ICOHP (eV/MO)" \
    -o ICOHP_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_ICOHP.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "ICOHP (eV/M-O)" \
    -o norm_ICOHP_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_norm_ICOHP.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "ICOBI (/MO)" \
    -o ICOBI_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_ICOBI.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "ICOBI (/M-O)" \
    -o norm_ICOBI_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_norm_ICOBI.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Gross population (Loewdin)" \
    -o GP_L_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_GP_Loewdin_M.tsv
    python ~/bin/tools/tetra/tsv.py -r ${row} -x "Metal (MO)" -y "Madelung energy (Loewdin, eV/MO)" \
    -o norm_Madelung_L_${row} /Users/hailey/Desktop/3_V_bulk/*_*_*/${row}/energy_norm_Madelung_Loewdin.tsv
done

### section5-2

cd /Users/hailey/Desktop/3_V_bulk/figures
python ~/bin/tools/tetra/concat.py -o energy --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_energy.tsv
python ~/bin/tools/tetra/concat.py -o norm_energy --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_norm_energy.tsv
python ~/bin/tools/tetra/concat.py -o norm_formation --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_norm_formation.tsv
python ~/bin/tools/tetra/concat.py -o norm_cohesive --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_norm_cohesive.tsv
python ~/bin/tools/tetra/concat.py -o ICOHP --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_ICOHP.tsv
python ~/bin/tools/tetra/concat.py -o norm_ICOHP --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_norm_ICOHP.tsv
python ~/bin/tools/tetra/concat.py -o ICOBI --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_ICOBI.tsv
python ~/bin/tools/tetra/concat.py -o norm_ICOBI --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_norm_ICOBI.tsv
python ~/bin/tools/tetra/concat.py -o bond --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_bond.tsv
python ~/bin/tools/tetra/concat.py -o norm_volume --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_norm_volume.tsv
python ~/bin/tools/tetra/concat.py -o chg --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_chg.tsv
python ~/bin/tools/tetra/concat.py -o mag --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_mag_M.tsv
python ~/bin/tools/tetra/concat.py -o norm_MadelungL --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_norm_Madelung_L.tsv
python ~/bin/tools/tetra/concat.py -o GrossPopulationL --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_GP_L_M.tsv
python ~/bin/tools/tetra/concat.py -o redoxP --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_redoxP.tsv
python ~/bin/tools/tetra/concat.py -o redoxP_clean --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_redoxP_clean.tsv
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
mv concat*rel.tsv rel6/
mv concat*.tsv rel7/

### section5-3

# cd /Users/hailey/Desktop/3_V_bulk/figures/rel6
# python ~/bin/tools/tetra/mendeleev2tsv.py -n 6 -p \
# group_id atomic_number atomic_volume  \
# boiling_point melting_point \
# mass density dipole_polarizability en_pauling \
# covalent_radius metallic_radius vdw_radius \
# evaporation_heat fusion_heat  heat_of_formation \
# ionenergies[1] ionenergies[2] ionenergies[3]
# python ~/bin/tools/tetra/operator.py -o + \
# -x concat_ionenergies_1.tsv \
# -y concat_ionenergies_2.tsv \
# -z concat_ionenergies_12.tsv
# python ~/bin/tools/tetra/operator.py -o + \
# -x concat_evaporation_heat.tsv \
# -y concat_fusion_heat.tsv \
# -z concat_sublimation_heat.tsv
# keys=("coord" "redoxP" "redoxP_clean")
# values=("merged_coord.tsv" "merged_redoxP.tsv" "merged_redoxP_clean.tsv")
# for i in "${!keys[@]}"; do
#     key="${keys[$i]}"
#     file="${values[$i]}"
#     python ~/bin/tools/tetra/concat.py -o "$key" --X \
#     /Users/hailey/Desktop/3_V_bulk/1_Tetrahedral_WZ/"$file" \
#     /Users/hailey/Desktop/3_V_bulk/2_Tetrahedral_ZB/"$file" \
#     /Users/hailey/Desktop/3_V_bulk/3_SquarePlanar_TN/"$file" \
#     /Users/hailey/Desktop/3_V_bulk/4_SquarePlanar_PD/"$file" \
#     /Users/hailey/Desktop/3_V_bulk/5_SquarePlanar_NB/"$file" \
#     /Users/hailey/Desktop/3_V_bulk/6_Octahedral_RS/"$file" \
#     /Users/hailey/Desktop/3_V_bulk/7_Pyramidal_LT/"$file"
# done
# keys=("element" "row")
# values=("merged_element.tsv" "merged_row.tsv")
# for i in "${!keys[@]}"; do
#     key="${keys[$i]}"
#     file="${values[$i]}"
#     python ~/bin/tools/tetra/concat.py -o "$key" --X \
#         /Users/hailey/Desktop/3_V_bulk/metal/"$file" \
#         /Users/hailey/Desktop/3_V_bulk/metal/"$file" \
#         /Users/hailey/Desktop/3_V_bulk/metal/"$file" \
#         /Users/hailey/Desktop/3_V_bulk/metal/"$file" \
#         /Users/hailey/Desktop/3_V_bulk/metal/"$file" \
#         /Users/hailey/Desktop/3_V_bulk/metal/"$file" \
#         /Users/hailey/Desktop/3_V_bulk/metal/"$file"
# done
for file in ~/bin/tools/tetra/png_rel/lr*.sh; do
    sh $file
done

# ## section5-4

# cd /Users/hailey/Desktop/3_V_bulk/figures/rel7
# python ~/bin/tools/tetra/mendeleev2tsv.py -n 7 -p \
# group_id atomic_number atomic_volume  \
# boiling_point melting_point \
# mass density dipole_polarizability en_pauling \
# covalent_radius metallic_radius vdw_radius \
# evaporation_heat fusion_heat  heat_of_formation \
# ionenergies[1] ionenergies[2] ionenergies[3]
# python ~/bin/tools/tetra/operator.py -o + \
# -x concat_ionenergies_1.tsv \
# -y concat_ionenergies_2.tsv \
# -z concat_ionenergies_12.tsv
# python ~/bin/tools/tetra/operator.py -o + \
# -x concat_evaporation_heat.tsv \
# -y concat_fusion_heat.tsv \
# -z concat_sublimation_heat.tsv
# python ~/bin/tools/tetra/concat.py -o coord --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_coord.tsv
# python ~/bin/tools/tetra/concat.py -o element --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_element.tsv
# python ~/bin/tools/tetra/concat.py -o row --X /Users/hailey/Desktop/3_V_bulk/nersc/*_*_*/merged_row.tsv
for file in ~/bin/tools/tetra/png/lr*.sh; do
    sh $file
    sh ~/bin/tools/tetra/png/lr-Ef.sh
done
# for file in ~/bin/tools/tetra/png/gbr*.sh; do
#     sh $file
# done
# for file in ~/bin/tools/tetra/png/nn*.sh; do
#     sh $file
# done