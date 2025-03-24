# Update the output directory to the desired location
output_dir = os.path.join(base_path, "figures")
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Re-run the data processing and saving
for label, (ads_folder, ref_folder) in folders.items():
    for spin_ads in spin_states:
        for spin_ref in spin_states:
            ads_data = read_energy_data(ads_folder, spin_ads)
            ref_data = read_energy_data(ref_folder, spin_ref)

            if ads_data is not None and ref_data is not None:
                merged = pd.merge(ads_data, ref_data, on="U", suffixes=("_ads", "_ref"))
                merged["ads_energy"] = merged["G_ads"] - merged["G_ref"]

                csv_filename = f"{label}_{spin_ref}_to_{spin_ads}.csv"
                csv_path = os.path.join(output_dir, csv_filename)
                merged[["U", "ads_energy"]].to_csv(csv_path, index=False)

                # Fitting and plotting
                try:
                    popt, _ = curve_fit(quadratic, merged["U"], merged["ads_energy"])
                    U_fit = np.linspace(merged["U"].min(), merged["U"].max(), 100)
                    energy_fit = quadratic(U_fit, *popt)

                    plt.figure()
                    plt.plot(merged["U"], merged["ads_energy"], 'o', label="Data")
                    plt.plot(U_fit, energy_fit, '-', label="Fit")
                    plt.xlabel("Applied Potential (V)")
                    plt.ylabel("Adsorption Energy (eV)")
                    plt.title(f"{label} ({spin_ref} → {spin_ads})")
                    plt.legend()
                    plt.tight_layout()

                    png_filename = f"{label}_{spin_ref}_to_{spin_ads}.png"
                    plt.savefig(os.path.join(output_dir, png_filename))
                    plt.close()
                except Exception as e:
                    print(f"Fitting failed for {label} {spin_ref} -> {spin_ads}: {e}")
