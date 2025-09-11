import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Scaling relationship for OER
a1, b1 = 1.63, 0.952  # OH → O
a2, b2 = 0.84, 3.14  # OH → OOH
a3, b3 = 1.0, 2.9  # OH → OOH

# Gas phase energies
h2 = -6.77149190
h2o = -14.23091949

zpeh2o, cvh2o, tsh2o = 0.558, 0.103, 0.675
zpeh2, cvh2, tsh2 = 0.268, 0.0905, 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2
go = gh2o - gh2 - 1.229 * 2

# Adsorbate corrections
zpeoh, cvoh, tsoh = 0.376, 0.042, 0.066
zpeo, cvo, tso = 0.064, 0.034, 0.060
zpeooh, cvooh, tsooh = 0.471, 0.077, 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh

def main():
    """Initialize energy data and plot the OER volcano diagram."""
    df = pd.DataFrame(columns=["ΔG_OH", "ΔG_O"])
    
    # Get data from O-OH_vs_Ueff.py calculations
    ueff_values, O_series, OH_series = get_o_oh_data()
    O_clean, O_Ocovered, O_from_VO = O_series
    OH_clean, OH_Ocovered, OH_from_VO = OH_series
    
    # Add data points for each Ueff value and configuration
    for i, ueff in enumerate(ueff_values):
        # Clean surface
        doh_clean = OH_clean[i]
        do_clean = O_clean[i]
        # df.loc[f"Ueff_{ueff}_clean", ["ΔG_OH", "ΔG_O"]] = [doh_clean, do_clean]
        
        # O-covered surface
        doh_covered = OH_Ocovered[i]
        do_covered = O_Ocovered[i]
        df.loc[f"Ueff_{ueff}_O_covered", ["ΔG_OH", "ΔG_O"]] = [doh_covered, do_covered]
        
        # From VO
        doh_vo = OH_from_VO[i]
        do_vo = O_from_VO[i]
        df.loc[f"Ueff_{ueff}_from_VO", ["ΔG_OH", "ΔG_O"]] = [doh_vo, do_vo]

    print(df)
    volcano_oer(df)

def get_o_oh_data():
    """Get O and OH adsorption energies from O-OH_vs_Ueff.py calculations."""
    from ase.io import read
    from pathlib import Path
    
    root = "~/Desktop/3_RuO2/1_RuO2_110"
    root_path = Path(root).expanduser()
    
    # Paths
    v_v_path = root_path / "1_V_V"
    v_oh_path = root_path / "2_V_OH"
    v_o_path = root_path / "3_V_O"
    o_ohcus_path = root_path / "4_O_OH"
    o_o_path = root_path / "5_O_O"
    vo_o_path = root_path / "6_VO_O"
    ho_o_path = root_path / "7_HO_O"
    
    ueff_values = [0, 1, 2, 3, 4]
    
    def get_energy_from_json(json_path):
        try:
            atoms = read(json_path)
            energy = atoms.get_total_energy()
            return energy
        except Exception as e:
            print(f"Error reading {json_path}: {e}")
            return None
    
    def load_series(base_path, label):
        series = []
        for u in ueff_values:
            json_path = base_path / f"{u}_" / "final_with_calculator.json"
            e = get_energy_from_json(json_path)
            if e is None:
                print(f"Skip {label} Ueff={u} (missing/failed)")
                return None
            series.append(e)
        return series
    
    energies_v_v = load_series(v_v_path, "V_V")
    energies_v_oh = load_series(v_oh_path, "V_OH")
    energies_v_o = load_series(v_o_path, "V_O")
    energies_o_ohcus = load_series(o_ohcus_path, "O_OHcus")
    energies_o_o = load_series(o_o_path, "O_O")
    energies_vo_o = load_series(vo_o_path, "VO_O")
    energies_ho_o = load_series(ho_o_path, "HO_O")
    
    if any(x is None for x in [energies_v_v, energies_v_oh, energies_v_o, energies_o_ohcus, energies_o_o, energies_vo_o, energies_ho_o]):
        return None
    
    # Calculate adsorption free energies
    O_clean = []
    O_Ocovered = []
    O_from_VO = []
    OH_clean = []
    OH_Ocovered = []
    OH_from_VO = []
    
    for i, ueff in enumerate(ueff_values):
        # Gas phase energies
        o_gas_energy = gh2o - gh2
        oh_gas_energy = gh2o - gh2/2
        
        # O adsorption energies
        O_clean.append(energies_v_o[i] + dgo - energies_v_v[i] - o_gas_energy)
        O_Ocovered.append(energies_o_o[i] + dgo - energies_v_o[i] - o_gas_energy)
        O_from_VO.append(energies_o_o[i] + dgo - energies_vo_o[i] - o_gas_energy)
        
        # OH adsorption energies
        OH_clean.append(energies_v_oh[i] + dgoh - energies_v_v[i] - oh_gas_energy)
        OH_Ocovered.append(energies_o_ohcus[i] + dgoh - energies_v_o[i] - oh_gas_energy)
        OH_from_VO.append(energies_ho_o[i] + dgoh - energies_vo_o[i] - oh_gas_energy)
    
    return (
        ueff_values,
        (O_clean, O_Ocovered, O_from_VO),
        (OH_clean, OH_Ocovered, OH_from_VO),
    )

def get_oer_activity1(doh, do):
    """Calculates OER activity based on adsorption free energies."""    

    # Use scaling relation for OOH
    dooh = a2 * doh + b2
    dooh_ = a3 * doh + b3

    # OER steps: * → *OH → *O → *OOH → * + O2
    dg1 = doh  # * → *OH
    dg2 = do - doh  # *OH → *O
    dg3 = dooh - do  # *O → *OOH
    dg4 = 4.92 - dooh  # *OOH → * + O2

    dg3_ = dooh_ - do
    dg4_ = 4.92 - dooh_
    
    # Overpotential is the maximum step
    op = max([dg2, dg3]) - 1.23
    op_ = max([dg2, dg3_]) - 1.23

    return op, op_

def get_oer_activity2(doh):
    """Calculates OER activity based on adsorption free energies."""    

    # Use scaling relation for OOH
    do = a1 * doh + b1
    dooh = a2 * doh + b2
    dooh_ = a3 * doh + b3

    # OER steps: * → *OH → *O → *OOH → * + O2
    dg1 = doh  # * → *OH
    dg2 = do - doh  # *OH → *O
    dg3 = dooh - do  # *O → *OOH
    dg4 = 4.92 - dooh  # *OOH → * + O2

    dg3_ = dooh_ - do
    dg4_ = 4.92 - dooh_
    
    # Overpotential is the maximum step
    op = max([dg2, dg3]) - 1.23
    op_ = max([dg2, dg3_]) - 1.23

    return op, op_

def get_oer_activity(x): # x = do - doh
    """Calculates OER activity based on adsorption free energies."""    

    doh = (x - b1) / (a1 - 1)
    do = a1 * doh + b1
    dooh = a2 * doh + b2
    dooh_ = a3 * doh + b3

    # OER steps: * → *OH → *O → *OOH → * + O2
    dg1 = doh  # * → *OH
    dg2 = do - doh  # *OH → *O
    dg3 = dooh - do  # *O → *OOH
    dg4 = 4.92 - dooh  # *OOH → * + O2

    dg3_ = dooh_ - do
    dg4_ = 4.92 - dooh_
    
    # Overpotential is the maximum step
    op = max([dg2, dg3]) - 1.23
    op_ = max([dg2, dg3_]) - 1.23

    return op, op_

def volcano_oer(df):
    """Plots the OER volcano plot."""
    xmin, xmax, xtick = 0.8, 2.1, 0.2
    ymin, ymax, ytick = -0.1, 0.9, 0.1
    
    xx = np.linspace(xmin-3, xmax+3, 100)

    # Calculate OER steps using scaling relations
    # doh = xx
    # do = a1 * xx + b1
    # dooh = a2 * xx + b2
    # dooh_ = a3 * xx + b3
    
    # OER steps: * → *OH → *O → *OOH → * + O2
    dg1 = xx  # * → *OH
    dg2 = (a1 - 1) * xx + b1  # *OH → *O
    x_intercept = (a1 - 1) * (1.23 - (b2 - b1)) / (a2 - a1) + b1
    dg3 = (a2 - a1) * xx + (b2 - b1)  # *O → *OOH
    dg3_ = (a3 - a1) * xx + (b3 - b1)
    dg4 = 4.92 - (a2 * xx + b2)  # *OOH → * + O2
    dg4_ = 4.92 - (a3 * xx + b3)

    plt.figure(figsize=(5.5, 5))

    plt.axhline(y=0.22, color='black', linestyle='-', alpha=0.7, linewidth=1)
    plt.axhline(y=0.31, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plt.text(0.9, 0.22-0.01, r'sRuO$_2$', fontsize=10, color='black')
    plt.text(0.9, 0.31-0.01, r'cRuO$_2$', fontsize=10, color='black')
    
    # Plot OER steps
    # plt.plot(dg2, dg1-1.23, linewidth=2.0, color='wheat', label=r"$\ast \rightarrow \ast \mathrm{OH}$")
    plt.plot(dg2, dg2-1.23, linewidth=2.0, color='tan', label=r"$\ast \mathrm{OH} \rightarrow \ast \mathrm{O}$")
    plt.plot(dg2, dg3-1.23, linewidth=2.0, color='lightsteelblue', label=r"$\ast \mathrm{O} \rightarrow \ast \mathrm{OOH}$")
    plt.plot(dg2, dg3_-1.23, linewidth=2.0, color='lavender', label=r"$\ast \mathrm{O} \rightarrow \ast \mathrm{OOH}$")
    plt.plot(dg2, dg4-1.23, linewidth=2.0, color='lightcoral', label=r"$\ast \mathrm{OOH} \rightarrow \ast + \mathrm{O}_2$")
    plt.plot(dg2, dg4_-1.23, linewidth=2.0, color='mistyrose', label=r"$\ast \mathrm{OOH} \rightarrow \ast + \mathrm{O}_2$")
    
    # # Plot limiting potential (maximum of all steps)
    # limiting_potential = np.maximum.reduce([dg1, dg2, dg3, dg4])-1.23
    # plt.plot(dg2, limiting_potential, linewidth=3.0, color='red', linestyle='--', alpha=0.8, label="Limiting Potential")
    
    # Scatter plot with labels
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    markers = ['s', '^']
    
    for i, txt in enumerate(df.index):
        # Use ΔG_O - ΔG_OH as x-axis for OER volcano plot
        x = df["ΔG_O"].iloc[i] - df["ΔG_OH"].iloc[i]
        y1, y2 = get_oer_activity(x)
        ueff_idx = i // 2  # Ueff index (0-4)
        config_idx = i % 2  # Configuration index (0=clean, 1=O_covered, 2=from_VO)
        
        color = colors[ueff_idx]
        marker = markers[config_idx]

        plt.scatter(x, y2, s=50, zorder=6, facecolor='white', edgecolor=color, marker=marker)
        plt.scatter(x, y1, s=50, zorder=6, color=color, marker=marker)
        
        # # Only label every few points to avoid clutter
        # if i % 2 == 0:  # Only label clean surface points
        #     plt.annotate(f'Ueff={ueff_idx}', (x, y1), xytext=(5, 0),
        #                  textcoords="offset points", ha="left", va='center', fontsize=8)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax) 
    plt.gca().invert_yaxis()  # y축을 역순으로 설정 (-10이 위, 10이 아래)
    plt.xticks(np.arange(xmin, xmax+xtick, xtick)) 
    plt.yticks(np.arange(ymin, ymax+ytick, ytick)) 
    plt.xlabel(r'$\Delta G_{O} - \Delta G_{OH}$ (eV)', fontsize=14)
    plt.ylabel("Overpotential (V)", fontsize=14)
    
    # 첫 번째 legend: OER steps
    legend1 = plt.legend(loc="lower right", bbox_to_anchor=(0.98, 0.02), fontsize=10, ncol=1, columnspacing=0.5, labelspacing=0.3, framealpha=1.0)
    plt.gca().add_artist(legend1)
    
    # 두 번째 legend: 데이터 포인트 (Ueff 값과 마커)
    from matplotlib.lines import Line2D
    legend_elements1 = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=f'Ueff={i}') for i, color in enumerate(colors)]
    legend_elements2 = [Line2D([0], [0], marker=marker, color='k', markerfacecolor='w', markersize=8, label=label) for marker, label in zip(markers, ['AEM', 'LOM'])]
    legend_elements = legend_elements1 + legend_elements2
    
    # Ueff legend
    legend2 = plt.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(0.98, 0.98), fontsize=9, ncol=1)
    plt.gca().add_artist(legend2)
    
    # Marker legend
    # legend3 = plt.legend(handles=legend_elements, loc="lower left", bbox_to_anchor=(0.02, 0.02), fontsize=9, ncol=1)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('volcano_oer1.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
