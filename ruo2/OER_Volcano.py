import os
import pandas as pd
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["axes.prop_cycle"] = plt.cycler(color=plt.get_cmap("tab20").colors)

# gas (전역 상수: 두 스크립트에서 동일 값 사용)
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.558
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2
gh = gh2/2
go = gh2o - gh2
goh = gh2o - gh
goo = 2*gh2o - 2*gh2
gooh = 2*gh2o - 3*gh

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

zpeooh = 0.471
cvooh = 0.077
tsooh = 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh
dgh = dgoh - dgo
dgoo = dgooh - dgh

def extract_energy_from_json(json_path):
    """
    JSON 파일에서 ASE Atoms 객체를 읽고 에너지를 추출하는 함수
    
    Args:
        json_path (str): JSON 파일 경로
        
    Returns:
        float: 에너지 값 (eV)
    """
    try:
        atoms = read(json_path)
        energy = atoms.get_total_energy()
        return energy
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None

def collect_energy_data(base_path):
    """
    모든 폴더와 서브폴더에서 에너지 데이터를 수집하는 함수
    
    Args:
        base_path (str): 기본 경로 (~/Desktop/3_RuO2/)
        
    Returns:
        dict: 정리된 에너지 데이터
    """
    # ~를 홈 디렉토리로 확장
    base_path = os.path.expanduser(base_path)
    folders = ["1_V_V", "2_V_OH", "3_O_V", "4_O_OH", "5_O_O", "6_O_OOH", "7_OO_OH"]
    semi_folders = ["0_V_V", "6_V_OH", "1_O_V", "2_O_OH", "3_O_O", "4_O_OOH", "5_OO_OH"]
    subfolders = ["0_", "1_", "2_", "3_", "4_"]
    energy_data = {}
    
    for i, folder in enumerate(folders):
        folder_path = os.path.join(base_path, "1_RuO2_Ueff", folder)
        if not os.path.exists(folder_path):
            print(f"폴더가 존재하지 않습니다: {folder_path}")
            continue
            
        energy_data[folder] = {}
        
        for subfolder in subfolders:
            subfolder_path = os.path.join(folder_path, subfolder)
            json_path = os.path.join(subfolder_path, "final_with_calculator.json")
            
            if os.path.exists(json_path):
                energy = extract_energy_from_json(json_path)
                energy_data[folder][subfolder] = energy
            else:
                print(f"JSON 파일이 존재하지 않습니다: {json_path}")
                energy_data[folder][subfolder] = None

        if folder == "5_O_O" or folder == "7_OO_OH":
            for ueff in range(1, 10):
                subfolder_path = os.path.join(folder_path, "between_2_and_3", f"{ueff}_")
                json_path = os.path.join(subfolder_path, "final_with_calculator.json")
                energy = extract_energy_from_json(json_path)
                energy_data[folder][f"2.{ueff}_"] = energy

        semi_folder = semi_folders[i]

        subfolder_path = os.path.join(base_path, "2_RuO2_OER", "7_RuO2", semi_folder)
        json_path = os.path.join(subfolder_path, "final_with_calculator.json")
        energy = extract_energy_from_json(json_path)
        energy_data[folder]["2.5_"] = energy
        energy_data[folder]["RuO2"] = energy

        re_forder_names = ["ReRuO2", "RuO2_1%", "RuO2_2%", "ReRuO2_1%", "ReRuO2_2%"]
        for f, re_folder in enumerate(["3_ReRuO2", "4_RuO2_strain_XRD", "5_RuO2_strain_ReO2", "6_ReRuO2_strain_1%", "7_ReRuO2_strain_2%"]):
            folder_name = re_forder_names[f]
            subfolder_path = os.path.join(base_path, "3_ReRuO2_OER", re_folder, semi_folder)
            json_path = os.path.join(subfolder_path, "final_with_calculator.json")
            energy = extract_energy_from_json(json_path)
            energy_data[folder][folder_name] = energy
        
    return energy_data

def create_dataframe(energy_data):
    """
    에너지 데이터를 DataFrame으로 변환하는 함수
    컬럼은 폴더명으로 하고 인덱스는 0, 1, 2, 3, 4로 설정
    
    Args:
        energy_data (dict): 에너지 데이터 딕셔너리
        
    Returns:
        pd.DataFrame: 정리된 데이터프레임 (전치된 형태)
    """
    # 데이터를 딕셔너리 형태로 변환
    data_dict = {}
    
    for folder, subfolder_data in energy_data.items():
        for subfolder, energy in subfolder_data.items():
            # 서브폴더명에서 숫자만 추출 (0_, 1_ -> 0, 1)
            row_name = subfolder.rstrip('_')
            if row_name not in data_dict:
                data_dict[row_name] = {}
            data_dict[row_name][folder] = energy
    
    # DataFrame 생성 (인덱스는 0, 1, 2, 3, 4, 컬럼은 폴더명)
    df = pd.DataFrame(data_dict).T
    
    # 인덱스 순서 정렬
    df = df.reindex(index=['0', '1', '2', '2.1', '2.2', '2.3', '2.4', '2.5', '2.6', '2.7', '2.8', '2.9', '3', '4', 'RuO2', 'ReRuO2', 'RuO2_1%', 'RuO2_2%', 'ReRuO2_1%', 'ReRuO2_2%'])
    # df = df.reindex(index=['0', '1', '2', '2.1', '2.2', '2.3', '2.4', '2.5', '2.6', '2.7', '2.8', '2.9', '3', '4', 'RuO2', 'ReRuO2', 'RuO2_1%', 'RuO2_2%'])
    
    return df

def main():
    """
    메인 실행 함수
    """
    # 기본 경로 설정
    base_path = "~/Desktop/3_RuO2"
    energy_data = collect_energy_data(base_path)
    df = create_dataframe(energy_data)

    print(df)

    # folders = ["1_V_V", "2_V_OH", "3_O_V", "4_O_OH", "5_O_O", "6_O_OOH", "7_OO_OH"]
    gibbs_correction_o_v = 0.058092 - 0.033706
    gibbs_correction_o_oh = 0.439910 - 0.128005
    gibbs_correction_o_o = 0.156862 - 0.100681
    gibbs_correction_oh_oo = 0.502349 - 0.144209
    gibbs_correction_o_ooh = 0.502349 - 0.144209
    gibbs_correction_v_oh = 0.364538 - 0.074708

    df["ΔG_O"] = df["3_O_V"] + gibbs_correction_o_v - df["1_V_V"] - go
    df["ΔG_OH"] = df["2_V_OH"] + gibbs_correction_v_oh - df["1_V_V"] - goh
    df["ΔG_O_O"] = df["5_O_O"] + gibbs_correction_o_o - df["3_O_V"] - gibbs_correction_o_v - go
    df["ΔG_O_OH"] = df["4_O_OH"] + gibbs_correction_o_oh - df["3_O_V"] - gibbs_correction_o_v - goh
    df["ΔG_O_OOH"] = df["ΔG_O_OH"] + 3.2
    # df["ΔG_O_OOH"] = df["6_O_OOH"] + gibbs_correction_o_ooh - df["3_O_V"] - gibbs_correction_o_v - goh
    # df["ΔG_OO_OH"] = df["7_OO_OH"] + gibbs_correction_oh_oo - df["3_O_V"] - gibbs_correction_o_v - goh
    
    df["ΔG1: *O→*O+*OH"] = df["4_O_OH"] + gibbs_correction_o_oh - goh - df["3_O_V"] - gibbs_correction_o_v
    df["ΔG2: *O+*OH→*O+*O"] = df["5_O_O"] + gibbs_correction_o_o - go - df["4_O_OH"] - gibbs_correction_o_oh + goh
    # df["ΔG3a: *O+*O→*O+*OOH"] = df["6_O_OOH"] + gibbs_correction_o_ooh - goh - df["5_O_O"] - gibbs_correction_o_o
    df["ΔG3a: *O+*O→*O+*OOH"] = df["ΔG_O_OOH"] - df["ΔG_O_O"]
    df["ΔG3d: *O+*O→*OO+*OH"] = df["7_OO_OH"] + gibbs_correction_oh_oo - goh - df["5_O_O"] - gibbs_correction_o_o
    df["ΔG4a: *O+*OOH→*O+O2"] = df["3_O_V"] + gibbs_correction_o_v - df["6_O_OOH"] - gibbs_correction_o_ooh + gooh + 4.92
    df["ΔG4d: *OO+*OH→*OH+O2"] = df["2_V_OH"] + gibbs_correction_v_oh - df["7_OO_OH"] - gibbs_correction_oh_oo + goo + 4.92
    df["ΔG5d: *OH→*O"] = df["3_O_V"] + gibbs_correction_o_v - go - df["2_V_OH"] - gibbs_correction_v_oh + goh
    # df["ΔG4"] = 4.92 - df["ΔG1"] - df["ΔG2"] - df["ΔG3"] - df["ΔG4"] - df["ΔG5"]

    print(df)

    # 'ΔG'가 포함되지 않은 컬럼 제거
    dg_columns = [c for c in df.columns if "ΔG" in c]
    # Ueff 관련 행만 선택 (RuO2, ReRuO2 등 제외)
    ueff_rows = [c for c in df.index if not any(x in c for x in ["RuO2", "ReRuO2"])]
    df_dg = df[dg_columns]
    df_scaling = df_dg.loc[ueff_rows]

    # Scaling relationship plot: x=(ΔG_O_O - ΔG_O_OH), y in [ΔG1..ΔG5]
    if "ΔG_O_O" in df_dg.columns and "ΔG_O_OH" in df_dg.columns:
        # x_series = (df["ΔG_O"] - df["ΔG_OH"]).rename("ΔG_O - ΔG_OH")
        x_series = (df_dg["ΔG_O_O"] - df_dg["ΔG_O_OH"]).rename("ΔG_O_O - ΔG_O_OH")
        y_columns = [c for c in ["ΔG1: *O→*O+*OH", "ΔG2: *O+*OH→*O+*O", "ΔG3a: *O+*O→*O+*OOH", "ΔG3d: *O+*O→*OO+*OH", "ΔG4a: *O+*OOH→*O+O2", "ΔG4d: *OO+*OH→*OH+O2", "ΔG5d: *OH→*O"] if c in df.columns]
        if len(y_columns) > 0:
            plt.figure(figsize=(6, 4))
            xmin, xmax = 1.0, 2.0
            ymin, ymax = 1.3, 1.9
            # xmin, xmax = 1.0, 2.0
            # ymin, ymax = 0.0, 2.0
            colors = ["C7", "C3", "C5", "C1", "C9", "C11", "C13"]
            for idx, ycol in enumerate(y_columns):
                xy = pd.concat([x_series, df_dg[ycol]], axis=1).dropna()
                if xy.empty:
                    continue
                xvals = xy[x_series.name].to_numpy()
                yvals = xy[ycol].to_numpy()
                # Ueff 관련 데이터는 흰색, RuO2/ReRuO2 관련 데이터는 색칠
                for i, (x, y) in enumerate(zip(xvals, yvals)):
                    row_name = xy.index[i]
                    if row_name == "RuO2" and idx == 3:
                        plt.scatter(x, y, color='black', zorder=10)
                    elif row_name == "ReRuO2" and idx == 3:
                        plt.scatter(x, y, color='red', zorder=10)
                    elif row_name == "RuO2_1%" and idx == 3:
                        plt.scatter(x, y, color='green', zorder=10)
                    elif row_name == "RuO2_2%" and idx == 3:
                        plt.scatter(x, y, color='orange', zorder=10)
                    # elif row_name == "ReRuO2_1%" and idx == 3:
                    #     plt.scatter(x, y, color='magenta', zorder=10)
                    # elif row_name == "ReRuO2_2%" and idx == 3:
                    #     plt.scatter(x, y, color='blue', zorder=10)
                    # else:
                    #     plt.scatter(x, y, facecolor='white', edgecolor=colors[idx % len(colors)], zorder=9)
                # print(idx,xvals, yvals)
                # if idx == 2 or idx == 4:
                #     plt.scatter(xvals[3], yvals[3], color=colors[idx % len(colors)], zorder=10)
                # else:
                #     plt.scatter(xvals[5], yvals[5], color=colors[idx % len(colors)], zorder=10)

                # # ΔG1에만 annotation 추가
                # if ycol == "ΔG1: *O→*O+*OH":
                #     for i, (x, y) in enumerate(zip(xvals, yvals)):
                #         plt.annotate(f'Ueff =\n{i} eV', (x, y), xytext=(0, -7), 
                #                 textcoords='offset points', fontsize=8, ha='center', va='top')

                # # ΔG3d에만 annotation 추가
                # if ycol == "ΔG3d: *O+*O→*OO+*OH":
                #     for i, (x, y) in enumerate(zip(xvals, yvals)):
                #         if i == 0 or i == 2 or i == 4:
                #             plt.annotate(f'{i} eV', (x, y), xytext=(0, 7), 
                #                     textcoords='offset points', fontsize=8, ha='center')
                #         elif i == 1:
                #             plt.annotate(f'Ueff = {i} eV', (x, y), xytext=(7, 0), 
                #                     textcoords='offset points', fontsize=8, ha='left', va='center')
                #         elif i == 3:
                #             plt.annotate(f'{i} eV', (x, y), xytext=(7, 0), 
                #                     textcoords='offset points', fontsize=8, ha='left', va='center')

                if xvals.size >= 2:
                    m, b = np.polyfit(xvals, yvals, 1)
                    xs = np.linspace(xmin, xmax, 100)
                    plt.plot(xs, m*xs + b, color=colors[idx % len(colors)], label=f"{ycol}")

            # plt.scatter(df_dg["ΔG2: *O+*OH→*O+*O"]['RuO2'], df_dg["ΔG2: *O+*OH→*O+*O"]['RuO2'], facecolor='white', edgecolor=colors[idx % len(colors)], zorder=10)

            plt.axhline(y=1.23+0.196, color='silver', linestyle='-', linewidth=1, zorder=0)
            plt.axhline(y=1.23+0.343, color='silver', linestyle='--', linewidth=1, zorder=0)
            plt.text(1.86, 1.23+0.196-0.012, r'Re-RuO$_2$', fontsize=10, color='silver')
            plt.text(1.86, 1.23+0.343-0.012, r'cRuO$_2$', fontsize=10, color='silver')
            plt.xlabel(r'$\Delta G_{O} - \Delta G_{OH}$ (eV)')
            plt.ylabel(r'$\Delta G_{max}$ (eV)')
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            plt.grid(True, alpha=0.1)
            plt.gca().invert_yaxis()
            # plt.legend(loc='lower right', bbox_to_anchor=(0.99, 0.01), fontsize=9, ncol=1)
            # plt.legend(loc='upper right', bbox_to_anchor=(0.99, 0.99), fontsize=9, ncol=1)
            # tight_layout 대신 수동으로 여백 조정
            plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.12)
            out_png = "OER_scaling.png"
            plt.savefig(out_png, dpi=300, bbox_inches='tight')
            plt.show()
            plt.close()
            print(f"스케일링 플롯 저장됨: {out_png}")

    print(df_dg)

    return df, energy_data

if __name__ == "__main__":
    df, energy_data = main()
