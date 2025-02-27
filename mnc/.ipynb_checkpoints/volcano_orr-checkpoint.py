import numpy as np # v1.26.4
import pandas as pd # v2.2.2
import matplotlib.pyplot as plt # v3.8.3

# Scaling relationship
a1, b1 = 1.79, 0.62  # OH → O
a2, b2 = 0.87, 3.08  # OH → OOH

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
    """Initialize energy data and plot the ORR volcano diagram."""
    df = pd.DataFrame(columns=["E_", "E_OH", "ΔE_OH", "ΔG_OH", "lp"])
    df.loc["haily_clean", ["E_", "E_OH"]] = [-280.17697237, -290.94763409999996]
    df.loc["haily_antipodalOH", ["E_", "E_OH"]] = [-290.94763409999996, -300.6264046]
    df.loc["haily_adjacentOH", ["E_", "E_OH"]] = [-290.94763409999996, -300.7382192]
    df.loc["haily_antipodalO", ["E_", "E_OH"]] = [-285.94477754, -295.52658608]
    df.loc["haily_adjacentO", ["E_", "E_OH"]] = [-285.94477754, -295.22826514]
    # df.loc["haily_antipodalOOH", ["E_", "E_OH"]] = [-295.18886945, -304.81150537999997]
    # df.loc["haily_adjacentOOH", ["E_", "E_OH"]] = [-295.18886945, -305.39491661]
    # df.loc["roman_SCAN(HS→LS)", ["ΔE_OH"]] = 0.51
    df.loc["roman_PBE+U4.3(HS→LS)", ["ΔE_OH"]] = 0.35
    # df.loc["roman_BELYP(HS→LS)", ["ΔE_OH"]] = 0.60
    df.loc["roman_HSE06(HS→LS)", ["ΔE_OH"]] = 0.53
    # df.loc["roman_PBE0(HS→LS)", ["ΔE_OH"]] = 0.52
    df.loc["roman_DMC(HS→LS)", ["ΔE_OH"]] = -0.11
    # df.loc["roman_DMC(HS→LS)", ["ΔG_OH"]] = -0.11

    # Iterate over rows to apply conditions
    for index, row in df.iterrows():
        if pd.notna(row["E_"]) and pd.notna(row["E_OH"]):
            df.loc[index, ["ΔG_OH", "lp"]] = get_orr_activity1(row["E_"], row["E_OH"])
        elif pd.notna(row["ΔE_OH"]):
            df.loc[index, ["ΔG_OH", "lp"]] = get_orr_activity2(row["ΔE_OH"])
        elif pd.notna(row["ΔG_OH"]):
            df.loc[index, "lp"] = get_orr_activity3(row["ΔG_OH"])

    print(df)
    # Plot volcano diagram
    volcano_orr(df)

def get_orr_activity1(clean, oh):
    """Calculates ORR activity based on adsorption energies."""
    doh = oh - clean - (gh2o - gh2 / 2) + dgoh
    do = a1 * doh + b1
    dooh = a2 * doh + b2
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    # dgmax = 1.23 - min(dg14)
    lp = min(dg14)
    return doh, lp
    
def get_orr_activity2(eoh):
    """Calculates ORR activity based on adsorption energies."""
    # eoh = oh - clean - (h2o - h2 / 2)
    # doh = oh - clean - ((h2o + zpeh2o - tsh2o + cvh2o) - (h2 + zpeh2 - tsh2 + cvh2) / 2) + dgoh
    doh = eoh - ((zpeh2o - tsh2o + cvh2o) - (zpeh2 - tsh2 + cvh2) / 2) + dgoh
    do = a1 * doh + b1
    dooh = a2 * doh + b2
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    # dgmax = 1.23 - min(dg14)
    lp = min(dg14)
    return doh, lp

def get_orr_activity3(doh):
    """Calculates ORR activity based on adsorption energies."""
    do = a1 * doh + b1
    dooh = a2 * doh + b2
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    # dgmax = 1.23 - min(dg14)
    lp = min(dg14)
    return lp

def volcano_orr(df):
    """Plots the ORR volcano plot."""
    xmin, xmax, xtick = -1.0, 3.0, 0.5
    ymin, ymax, ytick = -0.5, 1.5, 0.25
    
    xx = np.linspace(xmin, xmax, 100)

    doh = xx
    do = a1 * xx + b1
    dooh = a2 * xx + b2
    do2 = 4.92 - dooh
    
    dg1 = xx
    dg2 = (a1 - 1) * xx + b1
    dg3 = (a2 - a1) * xx + (b2 - b1)
    dg4 = 4.92 - (a2 * xx + b2)

    plt.figure(figsize=(5, 4), dpi=200)
    plt.plot(xx, dg4, linewidth=1.0, color='lightsteelblue', label=r"$\ast \rightarrow \ast \mathrm{OOH}$")
    plt.plot(xx, dg3, linewidth=1.0, color='lavender', label=r"$\ast \mathrm{OOH} \rightarrow \ast \mathrm{O}$")
    plt.plot(xx, dg2, linewidth=1.0, color='wheat', label=r"$\ast \mathrm{O} \rightarrow \ast \mathrm{OH}$")
    plt.plot(xx, dg1, linewidth=1.0, color='tan', label=r"$\ast \mathrm{OH} \rightarrow \ast$")
    
    # Scatter plot with labels
    for i, txt in enumerate(df.index):
        x, y = df["ΔG_OH"].iloc[i], df["lp"].iloc[i]
        plt.scatter(x, y, s=20, zorder=6)
        ha = "right" if x < 1 else "left"
        plt.annotate(txt, (x, y), xytext=(5 if x >= 1 else -5, 0),
                     textcoords="offset points", ha=ha, fontsize=8)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xticks(np.arange(xmin, xmax+xtick, xtick), fontsize=8) 
    plt.yticks(np.arange(ymin, ymax+ytick, ytick), fontsize=8) 
    plt.xlabel(r"$\mathrm{\Delta G_{OH}}$ (eV)")
    plt.ylabel("Limiting Potential (V)")
    plt.legend(loc="lower center", bbox_to_anchor=(0.5, 0.02), fontsize=8, ncol=2, columnspacing=0.5)
    plt.tight_layout()
    plt.savefig('volcano_orr.png')
    plt.show()

if __name__ == "__main__":
    main()
