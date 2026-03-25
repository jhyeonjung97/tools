root="/Users/jiuy97/Desktop/9_pourbaixGC/3_FeNC_gc"
ads_list=("H2Ogas_-0.5H2gas_star__OHstar" "H2Ogas_-1.0H2gas_star__Ostar" "2.0H2Ogas_-1.0H2gas_star__OHstarOHstar" "2.0H2Ogas_-2.0H2gas_star__OstarOstar" "2.0H2Ogas_-1.5H2gas_star__OHstarOstar")
for spin in HS IS LS; do
    spin_low=$(echo "${spin}" | tr '[:upper:]' '[:lower:]')
    for NELECT_delta in +1.0 +0.5 +0.0 -0.5 -1.0 -1.5 -2.0; do
        # clean
        NELECT_ref=$(grep "NELECT" "${root}/3_clean/${spin_low}/OUTCAR" | awk '{print $3}')
        NELECT=$(echo "scale=1; (${NELECT_ref} + ${NELECT_delta#+})/1" | bc)
        source_dirs="${root}/3_clean/${spin_low}/nelect_${NELECT}"
        destination="${root}/JungHybrid2026/VASP-6.4.2/PBE+U+VASPsol/FeN4C8-${spin}-ΔNELECT=${NELECT_delta}/001"
        mkdir -p "${destination}"
        cp "${source_dirs}/OUTCAR" "${destination}"
        cp "${source_dirs}/final_with_calculator.json" "${destination}/empty_slab.json"
        # OH
        NELECT_ref=$(grep "NELECT" "${root}/4_OH/${spin_low}/OUTCAR" | awk '{print $3}')
        NELECT=$(echo "scale=1; (${NELECT_ref} + ${NELECT_delta#+})/1" | bc)
        source_dirs="${root}/4_OH/${spin_low}/nelect_${NELECT}"
        destination="${root}/JungHybrid2026/VASP-6.4.2/PBE+U+VASPsol/FeN4C8-${spin}-ΔNELECT=${NELECT_delta}/001/H2Ogas_-0.5H2gas_star__OHstar"
        mkdir -p "${destination}"
        cp "${source_dirs}/OUTCAR" "${source_dirs}/final_with_calculator.json" "${destination}"
        # O
        NELECT_ref=$(grep "NELECT" "${root}/5_O/${spin_low}/OUTCAR" | awk '{print $3}')
        NELECT=$(echo "scale=1; (${NELECT_ref} + ${NELECT_delta#+})/1" | bc)
        source_dirs="${root}/5_O/${spin_low}/nelect_${NELECT}"
        destination="${root}/JungHybrid2026/VASP-6.4.2/PBE+U+VASPsol/FeN4C8-${spin}-ΔNELECT=${NELECT_delta}/001/H2Ogas_-1.0H2gas_star__Ostar"
        mkdir -p "${destination}"
        cp "${source_dirs}/OUTCAR" "${source_dirs}/final_with_calculator.json" "${destination}"
        # OH-OH
        NELECT_ref=$(grep "NELECT" "${root}/6_OH_OH/${spin_low}/OUTCAR" | awk '{print $3}')
        NELECT=$(echo "scale=1; (${NELECT_ref} + ${NELECT_delta#+})/1" | bc)
        source_dirs="${root}/6_OH_OH/${spin_low}/nelect_${NELECT}"
        destination="${root}/JungHybrid2026/VASP-6.4.2/PBE+U+VASPsol/FeN4C8-${spin}-ΔNELECT=${NELECT_delta}/001/2.0H2Ogas_-1.0H2gas_star__OHstarOHstar"
        mkdir -p "${destination}"
        cp "${source_dirs}/OUTCAR" "${source_dirs}/final_with_calculator.json" "${destination}"
        # OH-O
        NELECT_ref=$(grep "NELECT" "${root}/7_OH_O/${spin_low}/OUTCAR" | awk '{print $3}')
        NELECT=$(echo "scale=1; (${NELECT_ref} + ${NELECT_delta#+})/1" | bc)
        source_dirs="${root}/7_OH_O/${spin_low}/nelect_${NELECT}"
        destination="${root}/JungHybrid2026/VASP-6.4.2/PBE+U+VASPsol/FeN4C8-${spin}-ΔNELECT=${NELECT_delta}/001/2.0H2Ogas_-1.5H2gas_star__OHstarOstar"
        mkdir -p "${destination}"
        cp "${source_dirs}/OUTCAR" "${source_dirs}/final_with_calculator.json" "${destination}"
        # O-O
        NELECT_ref=$(grep "NELECT" "${root}/8_O_O/${spin_low}/OUTCAR" | awk '{print $3}')
        NELECT=$(echo "scale=1; (${NELECT_ref} + ${NELECT_delta#+})/1" | bc)
        source_dirs="${root}/8_O_O/${spin_low}/nelect_${NELECT}"
        destination="${root}/JungHybrid2026/VASP-6.4.2/PBE+U+VASPsol/FeN4C8-${spin}-ΔNELECT=${NELECT_delta}/001/2.0H2Ogas_-2.0H2gas_star__OstarOstar"
        mkdir -p "${destination}"
        cp "${source_dirs}/OUTCAR" "${source_dirs}/final_with_calculator.json" "${destination}"
    done
done