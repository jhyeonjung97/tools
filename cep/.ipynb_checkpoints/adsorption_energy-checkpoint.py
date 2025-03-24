import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Reference energies
h2 = -6.77149190
h2o = -14.23091949
o2 = 2 * h2o - 2 * h2 - 4.92
delta_G_OH = h2o - 0.5 * h2
delta_G_O = 0.5 * o2

# Folder settings
base_path = "/pscratch/sd/j/jiuy97/9_pourbaixGC"
output_dir = os.path.join(base_path, "figures")
os.makedirs(output_dir, exist_ok=True)

folders = {
    "OH_from_clean": ("4_OH", "3_clean"),
    "O_from_clean": ("5_O", "3_clean"),
    "OH_from_OH": ("6_OH_OH", "4_OH"),
    "O_from_OH": ("7_OH_O", "4_OH"),
    "OH_from_O": ("7_OH_O", "5_O"),
    "O_from_O": ("8_O_O", "5_O")
}

adsorption_refs = {
    "OH_from_clean": delta_G_OH,
    "OH_from_OH": delta_G_OH,
    "OH_from_O": delta_G_OH,
    "O_from_clean": delta_G_O,
    "O_from_OH": delta_G_O,
    "O_from_O": delta_G_O,
}

spin_states = ["hs", "is", "ls"]
spin_labels = {"hs": "HS", "is": "IS", "ls": "LS"}

# Plot settings
plot_order = [
    "HS→HS", "HS→IS", "HS→LS",
    "IS→HS", "IS→IS", "IS→LS",
    "LS→HS", "LS→IS", "LS→LS"
]
color_map = {
    "HS": ["#08306B", "#2171B5", "#DEEBF7"],
    "IS": ["#7F2704", "#FC9272", "#FEE0D2"],
    "LS": ["#00441B", "#41AB5D", "#C7E9C0"],
}
label_to_color = {}
for spin in ["HS", "IS", "LS"]:
    transitions = [f"{spin}→HS", f"{spin}→IS", f"{spin}→LS"]
    for i, label in enumerate(transitions):
        label_to_color[label] = color_map[spin][i]

# Fitting function
def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

# Read GCFE file
def read_energy_data(folder, spin):
    file_path = os.path.join(base_path, folder, spin, "GCFE_data_FULL.dat")
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["U", "G"])
        return df
    return None

# Main loop
for label, (ads_folder, ref_folder) in folders.items():
    ref_energy = adsorption_refs[label]
    all_data = []

    for spin_ads in spin_states:
        for spin_ref in spin_states:
            if (ads_folder in ["7_OH_O", "8_O_O"]) and spin_ads == "is":
                continue

            ads_data = read_energy_data(ads_folder, spin_ads)
            ref_data = read_energy_data(ref_folder, spin_ref)
            if ads_data is None or ref_data is None or ads_data.empty or ref_data.empty:
                continue

            try:
                popt_ads, _ = curve_fit(quadratic, ads_data["U"], ads_data["G"])
                popt_ref, _ = curve_fit(quadratic, ref_data["U"], ref_data["G"])
            except Exception as e:
                print(f"Fit failed for {label} {spin_ref}→{spin_ads}: {e}")
                continue

            U_common = np.linspace(-1.0, 2.0, 300)
            G_ads_fit = quadratic(U_common, *popt_ads)
            G_ref_fit = quadratic(U_common, *popt_ref)
            ΔG_ads = G_ads_fit - G_ref_fit - ref_energy

            label_name = f"{spin_labels[spin_ref]}→{spin_labels[spin_ads]}"
            for U, G in zip(U_common, ΔG_ads):
                all_data.append({"U": U, "corrected_ads_energy": G, "label": label_name})

    # Save combined CSV
    df_all = pd.DataFrame(all_data)
    df_all.to_csv(os.path.join(output_dir, f"{label}.csv"), index=False)

    # Plot all spin transitions in fixed order
    plt.figure()
    marker_U = [-0.5, 0.0, 0.5, 1.0, 1.5]

    for label_name in plot_order:
        df_subset = df_all[df_all["label"] == label_name]
        if df_subset.empty:
            continue
        try:
            popt, _ = curve_fit(quadratic, df_subset["U"], df_subset["corrected_ads_energy"])
            U_fit = np.linspace(-1.0, 2.0, 300)
            G_fit = quadratic(U_fit, *popt)
            color = label_to_color[label_name]
            plt.plot(U_fit, G_fit, label=label_name, color=color)

            for mu in marker_U:
                if -1.0 <= mu <= 2.0:
                    G_marker = quadratic(mu, *popt)
                    plt.scatter(mu, G_marker, color=color, edgecolors='k', zorder=5)
        except Exception as e:
            print(f"Final fit failed for {label} {label_name}: {e}")
            continue

    plt.xlabel("Applied Potential (V)")
    plt.ylabel("Corrected Adsorption Energy (eV)")
    plt.title(label.replace("_", " "))
    plt.xlim(-1.0, 2.0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{label}.png"))
    plt.close()
