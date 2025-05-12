import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def extract_mae_from_log(log_file):
    """로그 파일에서 test MAE 값을 추출하는 함수"""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            # test MAE: 0.xxxx 패턴을 찾음
            match = re.search(r'test MAE: (\d+\.\d+)', content)
            if match:
                return float(match.group(1))
    except FileNotFoundError:
        print(f"Warning: {log_file} not found")
        return None
    return None

def plot_mae_matrix():
    # 결과를 저장할 디렉토리
    root = '/Users/hailey/Desktop/7_V_bulk/figures'
    
    # ICOHP와 ionN의 정규화 방법
    norm_methods = ['x', 'n', 'o', 'c']  # x: no normalization, n: normalized, o: outer, c: core
    
    # 4x4 매트릭스를 저장할 배열
    mae_matrix = np.zeros((4, 4))
    
    # 각 조합에 대해 MAE 값을 추출
    for i, icohp_norm in enumerate(norm_methods):
        for j, ion_norm in enumerate(norm_methods):
            # 파일명 생성 (예: form_pred_cfse_gpr00_all_normxx.log)
            filename = f'form_pred_cfse_gpr00_all_norm{icohp_norm}{ion_norm}.log'
            filepath = os.path.join(root, filename)
            
            # MAE 값 추출
            mae = extract_mae_from_log(filepath)
            if mae is not None:
                mae_matrix[i, j] = mae
    
    # fig, ax를 명시적으로 생성 (figsize 6, 4)
    fig, ax = plt.subplots(figsize=(5, 3))
    sns.heatmap(
        mae_matrix,
        annot=True,
        fmt='.3f',
        cmap='Blues',
        xticklabels=norm_methods,
        yticklabels=norm_methods,
        cbar_kws={'label': 'Test MAE (eV)'},
        ax=ax
    )
    
    # xlabel, ylabel 변경
    plt.xlabel('Normalization for Ionization Energy')
    plt.ylabel('Normalization for -ICOHP')
    # 제목 제거 (아무것도 하지 않음)
    
    # x축과 y축 레이블에 설명 추가
    norm_labels = {
        'x': 'None',
        'n': 'OS*CN',
        'o': 'OS',
        'c': 'CN'
    }
    
    ax.set_xticks(np.arange(4) + 0.5)
    ax.set_xticklabels([norm_labels[m] for m in norm_methods])
    ax.set_yticks(np.arange(4) + 0.5)
    ax.set_yticklabels([norm_labels[m] for m in norm_methods], rotation=0)
    
    fig.tight_layout()
    
    # 결과 저장
    output_path = os.path.join(root, 'mae_matrix.png')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    print(f"MAE matrix plot saved as {output_path}")
    
    # 최적의 조합 찾기
    min_mae_idx = np.unravel_index(np.argmin(mae_matrix), mae_matrix.shape)
    best_icohp_norm = norm_methods[min_mae_idx[0]]
    best_ion_norm = norm_methods[min_mae_idx[1]]
    min_mae = mae_matrix[min_mae_idx]
    
    print(f"\n최적의 정규화 조합:")
    print(f"ICOHP: {norm_labels[best_icohp_norm]}")
    print(f"ionN: {norm_labels[best_ion_norm]}")
    print(f"Test MAE: {min_mae:.3f} eV")

if __name__ == '__main__':
    plot_mae_matrix() 