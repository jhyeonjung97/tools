import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define reference energies
h2 = -6.77149190
h2o = -14.23091949
ref_oh = h2o - h2/2
ref_o = h2o - h2 - 4.92/2

# Paths
base_path = "/pscratch/sd/j/jiuy97/9_pourbaixGC"
output_dir = os.path.join(base_path, "figures")
os.makedirs(output_dir, exist_ok=True)

# Folder and adsorption reference map
folders = {
    "OH_from_clean": ("4_OH", "3_clean"),
    "O_from_clean": ("5_O", "3_clean"),
    "OH_from_OH": ("6_OH_OH", "4_OH"),
    "O_from_OH": ("7_OH_O", "4_OH"),
    "OH_from_O": ("7_OH_O", "5_O"),
    "O_from_O": ("8_O_O", "5_O")
}
adsorption_refs = {
    "OH_from_clean": ref_oh,
    "OH_from_OH": ref_oh,
    "OH_from_O": ref_oh,
    "O_from_clean": ref_o,
    "O_from_OH": ref_o,
    "O_from_O": ref_o,
}
spin_states = ["hs", "is", "ls"]
spin_labels = {"hs": "HS", "is": "IS", "ls": "LS"}

# Fitting function
def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

# Read GCFE data
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

            U_min = max(ads_data["U"].min(), ref_data["U"].min())
            U_max = min(ads_data["U"].max(), ref_data["U"].max())
            if U_min >= U_max:
                print(f"Skipped {label} {spin_ref}→{spin_ads}: No overlapping U")
                continue

            U_common = np.linspace(U_min, U_max, 100)
            G_ads_fit = quadratic(U_common, *popt_ads)
            G_ref_fit = quadratic(U_common, *popt_ref)
            ΔG_ads = G_ads_fit - G_ref_fit - ref_energy

            spin_label = f"{spin_labels[spin_ref]}→{spin_labels[spin_ads]}"
            for U, G in zip(U_common, ΔG_ads):
                all_data.append({"U": U, "corrected_ads_energy": G, "label": spin_label})

    # Save combined CSV
    df_all = pd.DataFrame(all_data)
    df_all.to_csv(os.path.join(output_dir, f"{label}.csv"), index=False)

    # Plotting
    plt.figure()
    for spin_label in df_all["label"].unique():
        df_subset = df_all[df_all["label"] == spin_label]
        try:
            popt, _ = curve_fit(quadratic, df_subset["U"], df_subset["corrected_ads_energy"])
            U_fit = np.linspace(df_subset["U"].min(), df_subset["U"].max(), 200)
            G_fit = quadratic(U_fit, *popt)
            plt.plot(df_subset["U"], df_subset["corrected_ads_energy"], 'o', label=spin_label)
            plt.plot(U_fit, G_fit, '-', linewidth=1)
        except Exception as e:
            print(f"Final fit failed for {label} {spin_label}: {e}")
            continue

    plt.xlabel("Applied Potential (V)")
    plt.ylabel("Corrected Adsorption Energy (eV)")
    plt.title(label.replace("_", " "))
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{label}.png"))
    plt.close()
