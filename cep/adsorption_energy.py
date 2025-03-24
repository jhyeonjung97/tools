import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 기준 에너지 정의
h2 = -6.77149190
h2o = -14.23091949
o2 = 2 * h2o - 2 * h2 - 4.92
delta_G_OH = h2o - 0.5 * h2
delta_G_O = 0.5 * o2

# 경로 설정
base_path = "/pscratch/sd/j/jiuy97/9_pourbaixGC"
output_dir = os.path.join(base_path, "figures")
os.makedirs(output_dir, exist_ok=True)

# 폴더 매핑
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

def read_energy_data(folder, spin):
    file_path = os.path.join(base_path, folder, spin, "GCFE_data_FULL.dat")
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["U", "G"])
        print(df)
        return df
    else:
        return None

def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

for label, (ads_folder, ref_folder) in folders.items():
    for spin_ads in spin_states:
        for spin_ref in spin_states:
            if (ads_folder in ["7_OH_O", "8_O_O"]) and spin_ads == "is":
                continue

            ads_data = read_energy_data(ads_folder, spin_ads)
            ref_data = read_energy_data(ref_folder, spin_ref)

            if ads_data is not None and ref_data is not None:
                merged = pd.merge(ads_data, ref_data, on="U", suffixes=("_ads", "_ref"))
                merged["raw_ads_energy"] = merged["G_ads"] - merged["G_ref"]
                merged["corrected_ads_energy"] = merged["raw_ads_energy"] - adsorption_refs[label]

                csv_filename = f"{label}_{spin_ref}_to_{spin_ads}.csv"
                merged[["U", "corrected_ads_energy"]].to_csv(os.path.join(output_dir, csv_filename), index=False)

                try:
                    popt, _ = curve_fit(quadratic, merged["U"], merged["corrected_ads_energy"])
                    U_fit = np.linspace(merged["U"].min(), merged["U"].max(), 100)
                    energy_fit = quadratic(U_fit, *popt)

                    plt.figure()
                    plt.plot(merged["U"], merged["corrected_ads_energy"], 'o', label="Data")
                    plt.plot(U_fit, energy_fit, '-', label="Quadratic Fit")
                    plt.xlabel("Applied Potential (V)")
                    plt.ylabel("Corrected Adsorption Energy (eV)")
                    plt.title(f"{label} ({spin_ref} → {spin_ads})")
                    plt.legend()
                    plt.tight_layout()

                    png_filename = f"{label}_{spin_ref}_to_{spin_ads}.png"
                    plt.savefig(os.path.join(output_dir, png_filename))
                    plt.close()
                except Exception as e:
                    print(f"Fitting failed for {label} {spin_ref} -> {spin_ads}: {e}")
