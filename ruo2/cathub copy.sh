root="/Users/jiuy97/Desktop/3_RuO2"
cathub="${root}/JungHybrid2026/VASP-6.4.2/PBE+U+VASPsol"

for dir in ${root}/1_RuO2_Ueff/*_*/*_/; do
    ueff=$(basename "$dir" | cut -d'_' -f1)'.0'
    ads=$(basename "$(dirname "$dir")")
    if [[ "$ads" == "0_bulk" ]]; then
        mkdir -p "${cathub}/RuO2-Ueff=${ueff}"
        cp "${dir}/final_with_calculator.json" "${cathub}/RuO2-Ueff=${ueff}/bulk.json"
    elif [[ "$ads" == "1_V_V" ]]; then
        mkdir -p "${cathub}/RuO2-Ueff=${ueff}/110"
        cp "${dir}/final_with_calculator.json" "${cathub}/RuO2-Ueff=${ueff}/110/empty_slab.json"
    else
        if [[ "$ads" == "2_V_OH" ]]; then
            dirname="1.0H2Ogas_-0.5H2gas_star__OHstar"
        elif [[ "$ads" == "3_O_V" ]]; then
            dirname="1.0H2Ogas_-1.0H2gas_star__Ostar"
        elif [[ "$ads" == "4_O_OH" ]]; then
            dirname="2.0H2Ogas_-1.5H2gas_star__OstarOHstar"
        elif [[ "$ads" == "5_O_O" ]]; then
            dirname="2.0H2Ogas_-2.0H2gas_star__OstarOstar"
        else
            dirname="empty"
        fi
        if [[ "$dirname" != "empty" ]]; then
            mkdir -p "${cathub}/RuO2-Ueff=${ueff}/110/${dirname}"
            cp "${dir}/final_with_calculator.json" "${cathub}/RuO2-Ueff=${ueff}/110/${dirname}"
        fi
    fi 
done

mkdir -p "${cathub}/RuO2"
cp "${root}/1_RuO2_Ueff/0_bulk/between_2_and_3/5_/final_with_calculator.json" "${cathub}/RuO2/bulk.json"

for dir in ${root}/2_RuO2_OER/7_RuO2/*_*/; do
    ads=$(basename "$dir")
    if [[ "$ads" == "0_V_V" ]]; then
        mkdir -p "${cathub}/RuO2/110"
        cp "${dir}/final_with_calculator.json" "${cathub}/RuO2/110/empty_slab.json"
    else
        if [[ "$ads" == "1_O_V" ]]; then
            dirname="1.0H2Ogas_-1.0H2gas_star__Ostar"
        elif [[ "$ads" == "2_O_OH" ]]; then
            dirname="2.0H2Ogas_-1.5H2gas_star__OstarOHstar"
        elif [[ "$ads" == "3_O_O" ]]; then
            dirname="2.0H2Ogas_-2.0H2gas_star__OstarOstar"
        elif [[ "$ads" == "4_O_OOH" ]]; then
            dirname="3.0H2Ogas_-2.5H2gas_star__OstarOOHstar"
        elif [[ "$ads" == "5_OO_OH" ]]; then
            dirname="3.0H2Ogas_-2.5H2gas_star__OOstarOHstar"
        elif [[ "$ads" == "6_V_OH" ]]; then
            dirname="1.0H2Ogas_-0.5H2gas_star__OHstar"
        elif [[ "$ads" == "7_OH_OH" ]]; then
            dirname="2.0H2Ogas_-1.0H2gas_star__OHstarOHstar"
        else
            dirname="empty"
        fi
        if [[ "$dirname" != "empty" ]]; then
            mkdir -p "${cathub}/RuO2/110/${dirname}"
            cp "${dir}/final_with_calculator.json" "${cathub}/RuO2/110/${dirname}"
        fi
    fi 
done

for dir in ${root}/2_RuO2_OER/8_RuO2_brg/*_*/; do
    ads=$(basename "$dir")
    if [[ "$ads" == "1_H_V_V" ]]; then
        dirname="1.0H2gas_star__HstarHstar"
    elif [[ "$ads" == "2_H_O_O" ]]; then
        dirname="2.0H2Ogas_-1.0H2gas_star__OstarOstarHstarHstar"
    elif [[ "$ads" == "3_H_OH_OH" ]]; then
        dirname="2.0H2Ogas_star__OHstarOHstarHstarHstar"
    elif [[ "$ads" == "4_H_V_OH" ]]; then
        dirname="1.0H2Ogas_0.5H2gas_star__OHstarHstarHstar"
    elif [[ "$ads" == "5_H_V_O" ]]; then
        dirname="1.0H2Ogas_star__OstarHstarHstar"
    elif [[ "$ads" == "6_H_O_OH" ]]; then
        dirname="2.0H2Ogas_-0.5H2gas_star__OstarOHstarHstarHstar"
    else
        dirname="empty"
    fi
    if [[ "$dirname" != "empty" ]]; then
        mkdir -p "${cathub}/RuO2/110/${dirname}"
        cp "${dir}/final_with_calculator.json" "${cathub}/RuO2/110/${dirname}"
    fi 
done

dir="${root}/2_RuO2_OER/7_RuO2/3_O_O/"
dirname="0.5Ru2gas_2.0H2Ogas_-2.0H2gas_star__RuO2star"
mkdir -p "${cathub}/RuO2-O-vacRuO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/RuO2-O-vacRuO2/110/${dirname}"

dir="${root}/2_RuO2_OER/7_RuO2/a_O_brg/"
dirname="0.5Ru2gas_1.0H2Ogas_-1.0H2gas_star__RuOstar"
mkdir -p "${cathub}/RuO2-O-vacRuO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/RuO2-O-vacRuO2/110/${dirname}"

dir="${root}/2_RuO2_OER/7_RuO2/b_O2_brg/"
dirname="0.5Ru2gas_star__Rustar"
mkdir -p "${cathub}/RuO2-O-vacRuO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/RuO2-O-vacRuO2/110/${dirname}"

dir="${root}/2_RuO2_OER/7_RuO2/c_RuO2_brg/"
cp "${dir}/final_with_calculator.json" "${cathub}/RuO2-O-vacRuO2/110/empty_slab.json"

for dir in ${root}/3_ReRuO2_OER/3_ReRuO2/*_*/; do
    ads=$(basename "$dir")
    if [[ "$ads" == "2_O_OH" ]]; then
        dirname="1.0H2Ogas_-0.5H2gas_star__OHstar"
    elif [[ "$ads" == "3_O_O" ]]; then
        dirname="1.0H2Ogas_-1.0H2gas_star__Ostar"
    elif [[ "$ads" == "4_O_OOH" ]]; then
        dirname="2.0H2Ogas_-1.5H2gas_star__OOHstar"
    elif [[ "$ads" == "5_OO_OH" ]]; then
        dirname="2.0H2Ogas_-1.5H2gas_star__OOstarHstar"
    elif [[ "$ads" == "6_V_OH" ]]; then
        dirname="0.5H2gas_star__Hstar"
    else
        dirname="empty"
    fi
    if [[ "$dirname" != "empty" ]]; then
        mkdir -p "${cathub}/RuO2-O/110/${dirname}"
        cp "${root}/2_RuO2_OER/7_RuO2/${ads}/final_with_calculator.json" "${cathub}/RuO2-O/110/${dirname}"
        mkdir -p "${cathub}/ReRuO2-O-brg/110/${dirname}"
        cp "${root}/3_ReRuO2_OER/3_ReRuO2/${ads}/final_with_calculator.json" "${cathub}/ReRuO2-O-brg/110/${dirname}"
        mkdir -p "${cathub}/RuO2-O-lattice+1%/110/${dirname}"
        cp "${root}/3_ReRuO2_OER/4_RuO2_1%/${ads}/final_with_calculator.json" "${cathub}/RuO2-O-lattice+1%/110/${dirname}"
        mkdir -p "${cathub}/RuO2-O-lattice+2%/110/${dirname}"
        cp "${root}/3_ReRuO2_OER/5_RuO2_2%/${ads}/final_with_calculator.json" "${cathub}/RuO2-O-lattice+2%/110/${dirname}"
        mkdir -p "${cathub}/ReRuO2-O-alloy/110/${dirname}"
        cp "${root}/5_ReRuO2_alloy/1_OER/${ads}/final_with_calculator.json" "${cathub}/ReRuO2-O-alloy/110/${dirname}"
    fi 
done

cp ${root}/2_RuO2_OER/7_RuO2/1_O_V/final_with_calculator.json ${cathub}/RuO2-O/110/empty_slab.json
cp ${root}/3_ReRuO2_OER/3_ReRuO2/1_O_V/final_with_calculator.json ${cathub}/ReRuO2-O-brg/110/empty_slab.json
cp ${root}/3_ReRuO2_OER/4_RuO2_1%/1_O_V/final_with_calculator.json ${cathub}/RuO2-O-lattice+1%/110/empty_slab.json
cp ${root}/3_ReRuO2_OER/5_RuO2_2%/1_O_V/final_with_calculator.json ${cathub}/RuO2-O-lattice+2%/110/empty_slab.json
cp ${root}/5_ReRuO2_alloy/1_OER/1_O_V/final_with_calculator.json ${cathub}/ReRuO2-O-alloy/110/empty_slab.json

dir="${root}/5_ReRuO2_alloy/1_OER/3_O_O/"
dirname="0.5Ru2gas_2.0H2Ogas_-2.0H2gas_star__RuO2star"
mkdir -p "${cathub}/ReRuO2-O-alloy-vacRuO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/ReRuO2-O-alloy-vacRuO2/110/${dirname}"
dirname="0.5Re2gas_2.0H2Ogas_-2.0H2gas_star__ReO2star"
mkdir -p "${cathub}/ReRuO2-O-alloy-vacReO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/ReRuO2-O-alloy-vacReO2/110/${dirname}"

dir="${root}/5_ReRuO2_alloy/2_dissolution/c_O_brg/"
dirname="0.5Ru2gas_1.0H2Ogas_-1.0H2gas_star__RuOstar"
mkdir -p "${cathub}/ReRuO2-O-alloy-vacRuO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/ReRuO2-O-alloy-vacRuO2/110/${dirname}"
dirname="0.5Re2gas_1.0H2Ogas_-1.0H2gas_star__ReOstar"
mkdir -p "${cathub}/ReRuO2-O-alloy-vacReO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/ReRuO2-O-alloy-vacReO2/110/${dirname}"

dir="${root}/5_ReRuO2_alloy/2_dissolution/d_O2_brg/"
dirname="0.5Ru2gas_star__Rustar"
mkdir -p "${cathub}/ReRuO2-O-alloy-vacRuO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/ReRuO2-O-alloy-vacRuO2/110/${dirname}"
dirname="0.5Re2gas_star__Restar"
mkdir -p "${cathub}/ReRuO2-O-alloy-vacReO2/110/${dirname}"
cp "${dir}/final_with_calculator.json" "${cathub}/ReRuO2-O-alloy-vacReO2/110/${dirname}"

cp ${root}/5_ReRuO2_alloy/2_dissolution/8_O_RuO2_brg/final_with_calculator.json ${cathub}/ReRuO2-O-alloy-vacRuO2/110/empty_slab.json
cp ${root}/5_ReRuO2_alloy/2_dissolution/b_O_ReO2_brg/final_with_calculator.json ${cathub}/ReRuO2-O-alloy-vacReO2/110/empty_slab.json