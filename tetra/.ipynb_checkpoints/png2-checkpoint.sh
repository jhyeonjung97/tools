#!/bin/bash
if [[ ${here} == 'slac' ]]; then
    /usr/bin/rsync -e 'ssh' --ignore-times --size-only -avlzp -K --max-size=50000m \
    	--include="*/" \
    	--include="*.png" \
    	--exclude="*" \
     jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/3_V_bulk/nersc .
elif [[ ${here} == 'nersc' ]]; then
    
    ### section1
    
    # for dir in /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/*/
    # do
    #     cd $dir
    #     python ~/bin/tools/tetra/energy.py --save -x "Metal (M)" -y "Total energy (eV/M)" -n m
    # done
    # cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal
    # python ~/bin/tools/tetra/tsv2.py -l  3d 4d 5d -x "Metal (MO)" -y "Total energy (eV)" -o norm_energy *d/energy_norm_energy.tsv2
    # cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/oxide/0_min
    # python ~/bin/tools/tetra/energy.py --save -x "Metal (MxOy)" -y "Total energy (eV/M)" -n m
    # cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/6_Octahedral_RS
    # python ~/bin/tools/tetra/mendeleev2tsv2.py -p \
    # group_id atomic_number atomic_volume  \
    # boiling_point melting_point \
    # mass density dipole_polarizability en_pauling \
    # covalent_radius metallic_radius vdw_radius \
    # evaporation_heat fusion_heat heat_of_formation \
    # ionenergies[1] ionenergies[2] ionenergies[3]
    # python ~/bin/tools/tetra/operator.py -o + \
    # -x mendeleev_ionenergies_1.tsv2 \
    # -y mendeleev_ionenergies_2.tsv2 \
    # -z mendeleev_ionenergies_12.tsv2
    # python ~/bin/tools/tetra/operator.py -o + \
    # -x mendeleev_evaporation_heat.tsv2 \
    # -y mendeleev_fusion_heat.tsv2 \
    # -z mendeleev_sublimation_heat.tsv2
    # sed -i 's/\x0//g' *.tsv2
    
    ### section2
    
    # for dir in /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/*/
    # do
    #     cd $dir
    #     if [[ $dir == *'Tetrahedral'* ]]; then
    #         n=4; #python ~/bin/tools/tetra/energy.py --save -p hexa -x "Metal (MO)" -y "Hexagonal ratio [c/a]"
    #     elif [[ $dir == *'Pyramidal'* ]] || [[ $dir == *'Tetrahedral'* ]] || [[ $dir == *'SquarePlanar'* ]]; then
    #         n=4; #python ~/bin/tools/tetra/energy.py --save -p hexa -x "Metal (MO)" -y "Square prism ratio [c/a]"
    #     elif [[ $dir == *'Octahedral'* ]]; then
    #         n=6
    #     fi
    #     python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV)"
    #     python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV/MO)" -n m
    #     python ~/bin/tools/tetra/energy.py --save -p bond  -x "Metal (MO)" -y "Bond length (A)"
    #     python ~/bin/tools/tetra/energy.py --save -p bond -x "Metal (MO)" -y "Bond length (A/M-O)" -n $n
    #     python ~/bin/tools/tetra/energy.py --save -p volume -x "Metal (MO)" -y "Volume (A^3/MO)" -n m
    #     python ~/bin/tools/tetra/energy.py --save -p chg -e M  -x "Metal (MO)" -y "Bader charge (e-)"
    #     python ~/bin/tools/tetra/energy.py --save -p mag -e M -x "Metal (MO)" -y "|Magnetization|"
    #     python ~/bin/tools/tetra/energy.py --save -p ICOHP -x "Metal (MO)" -y "ICOHP (eV/MO)"
    #     python ~/bin/tools/tetra/energy.py --save -p ICOHP -x "Metal (MO)" -y "ICOHP (eV/M-O)" -n $n
    #     python ~/bin/tools/tetra/energy.py --save -p ICOBI -x "Metal (MO)" -y "ICOBI (/MO)"
    #     python ~/bin/tools/tetra/energy.py --save -p ICOBI -x "Metal (MO)" -y "ICOBI (eV/M-O)" -n $n
    #     python ~/bin/tools/tetra/energy.py --save -p GP_L -e M  -x "Metal (MO)" -y "Gross population (Loewdin)"
    #     python ~/bin/tools/tetra/energy.py --save -p Madelung_L -x "Metal (MO)" -y "Madelugn energy (Loewdin, eV/MO)" -n m
    #     python ~/bin/tools/tetra/formation_energy.py
    #     python ~/bin/tools/tetra/cohesive_energy.py
    #     # python ~/bin/tools/tetra/energy.py --save -p PSCENC -x "Metal (MO)" -y "PSCENC (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p TEWEN -x "Metal (MO)" -y "TEWEN (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p DENC -x "Metal (MO)" -y "DENC (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EXHF -x "Metal (MO)" -y "EXHF (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p XCENC -x "Metal (MO)" -y "XCENC (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p PAW_double_counting -x "Metal (MO)" -y "PAW_double_counting (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EENTRO -x "Metal (MO)" -y "EENTRO (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EBANDS -x "Metal (MO)" -y "EBANDS (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EATOM -x "Metal (MO)" -y "EATOM (eV/MO)" -n m
    #     sed -i 's/\x0//g' *.tsv2
    # done
    
    ### section3
    
    # for dir in /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/; do
    #     cd $dir
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Total energy (eV)" \
    #     -o energy 3d/energy_energy.tsv2 4d/energy_energy.tsv2 5d/energy_energy.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Total energy (eV/MO)" \
    #     -o norm_energy 3d/energy_norm_energy.tsv2 4d/energy_norm_energy.tsv2 5d/energy_norm_energy.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Bond length (A/M-O)" \
    #     -o bond 3d/energy_bond.tsv2 4d/energy_bond.tsv2 5d/energy_bond.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Volume (A^3/MO)" \
    #     -o norm_volume 3d/energy_norm_volume.tsv2 4d/energy_norm_volume.tsv2 5d/energy_norm_volume.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Bader charge (e-)" \
    #     -o chg 3d/energy_chg_M.tsv2 4d/energy_chg_M.tsv2 5d/energy_chg_M.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "|Magnetization|" \
    #     -o mag_M 3d/energy_mag_M.tsv2 4d/energy_mag_M.tsv2 5d/energy_mag_M.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOHP (eV/MO)" \
    #     -o ICOHP 3d/energy_ICOHP.tsv2 4d/energy_ICOHP.tsv2 5d/energy_ICOHP.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOHP (eV/M-O)" \
    #     -o norm_ICOHP 3d/energy_norm_ICOHP.tsv2 4d/energy_norm_ICOHP.tsv2 5d/energy_norm_ICOHP.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOBI (/MO)" \
    #     -o ICOBI 3d/energy_ICOBI.tsv2 4d/energy_ICOBI.tsv2 5d/energy_ICOBI.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "ICOBI (/M-O)" \
    #     -o norm_ICOBI 3d/energy_norm_ICOBI.tsv2 4d/energy_norm_ICOBI.tsv2 5d/energy_norm_ICOBI.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Gross population (Loewdin)" \
    #     -o GP_L_M 3d/energy_GP_Loewdin_M.tsv2 4d/energy_GP_Loewdin_M.tsv2 5d/energy_GP_Loewdin_M.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Madelung energy (Loewdin, eV/MO)" \
    #     -o norm_Madelung_L 3d/energy_norm_Madelung_Loewdin.tsv2 4d/energy_norm_Madelung_Loewdin.tsv2 5d/energy_norm_Madelung_Loewdin.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Formation energy (eV/MO)" \
    #     -o norm_formation 3d/energy_norm_formation.tsv2 4d/energy_norm_formation.tsv2 5d/energy_norm_formation.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Cohesive energy (eV/MO)" \
    #     -o norm_cohesive 3d/energy_norm_cohesive.tsv2 4d/energy_norm_cohesive.tsv2 5d/energy_norm_cohesive.tsv2
    #     python ~/bin/tools/tetra/tsv2.py -l 3d fm -x "Metal (MO)" -y "Formation energy (eV/MO)" \
    #     -o AFMvsFM 3d/energy_norm_formation.tsv2 fm/energy_norm_formation.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "PSCENC (eV/MO)" \
    #     # -o norm_PSCENC 3d/energy_norm_PSCENC.tsv2 4d/energy_norm_PSCENC.tsv2 5d/energy_norm_PSCENC.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "TEWEN (eV/MO)" \
    #     # -o norm_TEWEN 3d/energy_norm_TEWEN.tsv2 4d/energy_norm_TEWEN.tsv2 5d/energy_norm_TEWEN.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "DENC (eV/MO)" \
    #     # -o norm_DENC 3d/energy_norm_DENC.tsv2 4d/energy_norm_DENC.tsv2 5d/energy_norm_DENC.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "EXHF (eV/MO)" \
    #     # -o norm_EXHF 3d/energy_norm_EXHF.tsv2 4d/energy_norm_EXHF.tsv2 5d/energy_norm_EXHF.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "XCENC (eV/MO)" \
    #     # -o norm_XCENC 3d/energy_norm_XCENC.tsv2 4d/energy_norm_XCENC.tsv2 5d/energy_norm_XCENC.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "PAW_double_counting (eV/MO)" \
    #     # -o norm_PAW_double_counting 3d/energy_norm_PAW_double_counting.tsv2 4d/energy_norm_PAW_double_counting.tsv2 5d/energy_norm_PAW_double_counting.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "EENTRO (eV/MO)" \
    #     # -o norm_EENTRO 3d/energy_norm_EENTRO.tsv2 4d/energy_norm_EENTRO.tsv2 5d/energy_norm_EENTRO.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "EBANDS (eV/MO)" \
    #     # -o norm_EBANDS 3d/energy_norm_EBANDS.tsv2 4d/energy_norm_EBANDS.tsv2 5d/energy_norm_EBANDS.tsv2
    #     # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "EATOM (eV/MO)" \
    #     # -o norm_EATOM 3d/energy_norm_EATOM.tsv2 4d/energy_norm_EATOM.tsv2 5d/energy_norm_EATOM.tsv2
    # done
    
    ### section4-1
    
    # cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/6_Octahedral_RS
    # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Standard reduction potential (V)" -o redoxP \
    # ../oxide/energy_redoxP_3d.tsv2 ../oxide/energy_redoxP_4d.tsv2 ../oxide/energy_redoxP_5d.tsv2
    # python ~/bin/tools/tetra/tsv2.py -l 3d 4d 5d -x "Metal (MO)" -y "Standard reduction potential (V)" -o redoxP_clean \
    # ../oxide/energy_redoxP_clean_3d.tsv2 ../oxide/energy_redoxP_clean_4d.tsv2 ../oxide/energy_redoxP_clean_5d.tsv2

    ### section4-2
    
    # cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc
    # sh ~/bin/verve/spread.sh 6_Octahedral_RS/merged_redoxP.tsv2
    # sh ~/bin/verve/spread.sh 6_Octahedral_RS/merged_redoxP_clean.tsv2
    # sh ~/bin/verve/spread.sh oxide/merged_coord.tsv2
    # sh ~/bin/verve/spread.sh metal/merged_element.tsv2
    # sh ~/bin/verve/spread.sh metal/merged_row.tsv2
    # sed -i -e 's/CN/WZ/g' 1_Tetrahedral_WZ/merged_coord.tsv2
    # sed -i -e 's/CN/ZB/g' 2_Tetrahedral_ZB/merged_coord.tsv2
    # sed -i -e 's/CN/TN/g' 3_Pyramidal_LT/merged_coord.tsv2
    # sed -i -e 's/CN/PD/g' 4_Square_Planar_TN/merged_coord.tsv2
    # sed -i -e 's/CN/NB/g' 5_Square_Planar_NB/merged_coord.tsv2
    # sed -i -e 's/CN/RS/g' 6_Octahedral_RS/merged_coord.tsv2
    
    ### section5-1
    
    cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/figures
    for row in 3d 4d 5d
    do
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Total energy (eV)" \
        -o energy_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_energy.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Total energy (eV/MO)" \
        -o norm_energy_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_norm_energy.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Formation energy (eV/MO)" \
        -o norm_formation_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_norm_formation.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Cohesive energy (eV/MO)" \
        -o norm_cohesive_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_norm_cohesive.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Bond length (A/M-O)" \
        -o bond_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_bond.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Volume (A^3/MO)" \
        -o norm_volume_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_norm_volume.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Bader charge (e-)" \
        -o chg_M_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_chg_M.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "|Magnetization|" \
        -o mag_M_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_mag_M.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "ICOHP (eV/MO)" \
        -o ICOHP_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_ICOHP.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "ICOHP (eV/M-O)" \
        -o norm_ICOHP_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_norm_ICOHP.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "ICOBI (/MO)" \
        -o ICOBI_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_ICOBI.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "ICOBI (/M-O)" \
        -o norm_ICOBI_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_norm_ICOBI.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Gross population (Loewdin)" \
        -o GP_L_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_GP_Loewdin_M.tsv2
        python ~/bin/tools/tetra/tsv2.py -r ${row} -x "Metal (MO)" -y "Madelung energy (Loewdin, eV/MO)" \
        -o norm_Madelung_L_${row} /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/${row}/energy_norm_Madelung_Loewdin.tsv2
    done
    
    ### section5-2
    
    cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/figures
    python ~/bin/tools/tetra/concat.py -o energy --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_energy.tsv2
    python ~/bin/tools/tetra/concat.py -o norm_energy --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_norm_energy.tsv2
    python ~/bin/tools/tetra/concat.py -o norm_formation --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_norm_formation.tsv2
    python ~/bin/tools/tetra/concat.py -o norm_cohesive --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_norm_cohesive.tsv2
    python ~/bin/tools/tetra/concat.py -o ICOHP --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_ICOHP.tsv2
    python ~/bin/tools/tetra/concat.py -o norm_ICOHP --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_norm_ICOHP.tsv2
    python ~/bin/tools/tetra/concat.py -o ICOBI --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_ICOBI.tsv2
    python ~/bin/tools/tetra/concat.py -o norm_ICOBI --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_norm_ICOBI.tsv2
    python ~/bin/tools/tetra/concat.py -o bond --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_bond.tsv2
    python ~/bin/tools/tetra/concat.py -o norm_volume --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_norm_volume.tsv2
    python ~/bin/tools/tetra/concat.py -o chg --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_chg.tsv2
    python ~/bin/tools/tetra/concat.py -o mag --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_mag_M.tsv2
    python ~/bin/tools/tetra/concat.py -o norm_MadelungL --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_norm_Madelung_L.tsv2
    python ~/bin/tools/tetra/concat.py -o GrossPopulationL --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_GP_L_M.tsv2
    python ~/bin/tools/tetra/concat.py -o redoxP --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_redoxP.tsv2
    python ~/bin/tools/tetra/concat.py -o redoxP_clean --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_redoxP_clean.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_bond.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_chg.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_energy.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_GrossPopulationL.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_ICOBI.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_ICOHP.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_mag.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_norm_energy.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_norm_formation.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_norm_cohesive.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_norm_ICOHP.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_norm_ICOBI.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_norm_MadelungL.tsv2
    python ~/bin/tools/tetra/rel2octa.py concat_norm_volume.tsv2
    mv concat*rel.tsv2 rel8/
    mv concat*.tsv2 rel9/

    ### section5-3
    
    cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/figures/rel8
    python ~/bin/tools/tetra/mendeleev2tsv2.py -n 8 -p \
    group_id atomic_number atomic_volume  \
    boiling_point melting_point \
    mass density dipole_polarizability en_pauling \
    covalent_radius metallic_radius vdw_radius \
    evaporation_heat fusion_heat  heat_of_formation \
    ionenergies[1] ionenergies[2] ionenergies[3]
    python ~/bin/tools/tetra/operator.py -o + \
    -x concat_ionenergies_1.tsv2 \
    -y concat_ionenergies_2.tsv2 \
    -z concat_ionenergies_12.tsv2
    python ~/bin/tools/tetra/operator.py -o + \
    -x concat_evaporation_heat.tsv2 \
    -y concat_fusion_heat.tsv2 \
    -z concat_sublimation_heat.tsv2
    declare -A files_A
    files_A[coord]="merged_coord.tsv2"
    files_A[redoxP]="merged_redoxP.tsv2"
    files_A[redoxP_clean]="merged_redoxP_clean.tsv2"
    for key in "${!files_A[@]}"; do
        python ~/bin/tools/tetra/concat.py -o $key --X \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/1_Tetrahedral_WZ/${files_A[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/2_Tetrahedral_ZB/${files_A[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/3_SquarePlanar_TN/${files_A[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/4_SquarePlanar_PD/${files_A[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/5_SquarePlanar_NB/${files_A[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/7_Pyramidal_LT/${files_A[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/8_Tetrahedral_AQ/${files_A[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/9_SquarePlanar_AU/${files_A[$key]}
    done
    declare -A files_B
    files_B[element]="merged_element.tsv2"
    files_B[row]="merged_row.tsv2"
    for key in "${!files_B[@]}"; do
        python ~/bin/tools/tetra/concat.py -o $key --X \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]} \
        /pscratch/sd/j/jiuy97/3_V_bulk/nersc/metal/${files_B[$key]}
    done
    for file in ~/bin/tools/tetra/png_rel/lr*.sh; do
        sh $file
    done
    # for file in ~/bin/tools/tetra/png_rel/gbr*.sh; do
    #     sh $file
    # done
    # for file in ~/bin/tools/tetra/png_rel/nn*.sh; do
    #     sh $file
    # done
    
    ### section5-4
    
    cd /pscratch/sd/j/jiuy97/3_V_bulk/nersc/figures/rel9
    python ~/bin/tools/tetra/mendeleev2tsv2.py -n 9 -p \
    group_id atomic_number atomic_volume  \
    boiling_point melting_point \
    mass density dipole_polarizability en_pauling \
    covalent_radius metallic_radius vdw_radius \
    evaporation_heat fusion_heat  heat_of_formation \
    ionenergies[1] ionenergies[2] ionenergies[3]
    python ~/bin/tools/tetra/operator.py -o + \
    -x concat_ionenergies_1.tsv2 \
    -y concat_ionenergies_2.tsv2 \
    -z concat_ionenergies_12.tsv2
    python ~/bin/tools/tetra/operator.py -o + \
    -x concat_evaporation_heat.tsv2 \
    -y concat_fusion_heat.tsv2 \
    -z concat_sublimation_heat.tsv2
    python ~/bin/tools/tetra/concat.py -o coord --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_coord.tsv2
    python ~/bin/tools/tetra/concat.py -o element --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_element.tsv2
    python ~/bin/tools/tetra/concat.py -o row --X /pscratch/sd/j/jiuy97/3_V_bulk/nersc/*_*_*/merged_row.tsv2
    for file in ~/bin/tools/tetra/png/lr*.sh; do
        sh $file
    done
    # for file in ~/bin/tools/tetra/png/gbr*.sh; do
    #     sh $file
    # done
    # for file in ~/bin/tools/tetra/png/nn*.sh; do
    #     sh $file
    # done
    
    ### section6
    
    # for dir in /pscratch/sd/j/jiuy97/4_V_slab/*_*_*/*d/; do
    #     cd $dir
    #     python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV)"
    #     python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV/MO)" -n m
    #     python ~/bin/tools/tetra/energy.py --save -p area -x "Metal (MO)" -y "Area (A^2)"
    #     # python ~/bin/tools/tetra/energy.py --save -p PSCENC -x "Metal (MO)" -y "PSCENC (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p TEWEN -x "Metal (MO)" -y "TEWEN (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p DENC -x "Metal (MO)" -y "DENC (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EXHF -x "Metal (MO)" -y "EXHF (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p XCENC -x "Metal (MO)" -y "XCENC (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p PAW_double_counting -x "Metal (MO)" -y "PAW_double_counting (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EENTRO -x "Metal (MO)" -y "EENTRO (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EBANDS -x "Metal (MO)" -y "EBANDS (eV/MO)" -n m
    #     # python ~/bin/tools/tetra/energy.py --save -p EATOM -x "Metal (MO)" -y "EATOM (eV/MO)" -n m
    #     sed -i 's/\x0//g' *.tsv2
    # done
fi