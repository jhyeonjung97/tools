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

            if ads_data is not None and ref_data is not None and not ads_data.empty and not ref_data.empty:
                try:
                    # 각각 개별적으로 fitting
                    popt_ads, _ = curve_fit(quadratic, ads_data["U"], ads_data["G"])
                    popt_ref, _ = curve_fit(quadratic, ref_data["U"], ref_data["G"])
                except Exception as e:
                    print(f"Fit failed for {label} {spin_ref}->{spin_ads}: {e}")
                    continue

                # 공통 potential 범위
                U_min = max(ads_data["U"].min(), ref_data["U"].min())
                U_max = min(ads_data["U"].max(), ref_data["U"].max())
                if U_min >= U_max:
                    print(f"Skipped {label} {spin_ref}->{spin_ads}: No overlapping U range")
                    continue

                U_common = np.linspace(U_min, U_max, 100)
                G_ads_fit = quadratic(U_common, *popt_ads)
                G_ref_fit = quadratic(U_common, *popt_ref)

                ΔG_ads = G_ads_fit - G_ref_fit - adsorption_refs[label]
                df_result = pd.DataFrame({"U": U_common, "corrected_ads_energy": ΔG_ads})

                # 저장
                csv_name = f"{label}_{spin_ref}_to_{spin_ads}.csv"
                df_result.to_csv(os.path.join(output_dir, csv_name), index=False)

                try:
                    # ΔG_ads 도 다시 fitting
                    popt_final, _ = curve_fit(quadratic, df_result["U"], df_result["corrected_ads_energy"])
                    U_fit = np.linspace(U_min, U_max, 200)
                    fit_energy = quadratic(U_fit, *popt_final)

                    plt.figure()
                    plt.plot(df_result["U"], df_result["corrected_ads_energy"], 'o', label="ΔG_ads")
                    plt.plot(U_fit, fit_energy, '-', label="Fit")
                    plt.xlabel("Applied Potential (V)")
                    plt.ylabel("Corrected Adsorption Energy (eV)")
                    plt.title(f"{label}: {spin_ref} → {spin_ads}")
                    plt.legend()
                    plt.tight_layout()

                    png_name = f"{label}_{spin_ref}_to_{spin_ads}.png"
                    plt.savefig(os.path.join(output_dir, png_name))
                    plt.close()
                except Exception as e:
                    print(f"Final fit failed for {label} {spin_ref}->{spin_ads}: {e}")
