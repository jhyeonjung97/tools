#!/bin/env python3
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D

data = [
    # [catal, e_, e_oh, e_o, e_ooh, e_o2]
    ['Mn-N-C', -281.65, -292.61, -287.44, -296.92, -292.44], 
    ['Fe-N-C', -280.18, -290.95, -285.94, -295.19, -290.62],
    ['Co-N-C', -279.55, -289.62, -284.79, -294.19, -289.90],
    ['Ni-N-C', -277.21, -287.28, -281.31, -291.72, -287.37],
]

# gas
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
go2 = 2*gh2o - 2*gh2 + 4.92

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
dgo2 = dgooh - dgh

steps = ['O2->OOH*', 'OOH*->O*', 'O*->OH*', 'OH*->H2O']

# CSV 파일로 저장 (소숫점 모두 포함)
with open('orr_overpotential_results.csv', 'w', newline='') as csvfile:
    fieldnames = ['catal', 'goh', 'go', 'gooh', 'dg0', 'dg1', 'dg2', 'dg3', 'dg4', 'op14', 'op04']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')
    
    # 헤더 작성
    writer.writeheader()
    
    for row in data:
        catal, e_, e_oh, e_o, e_ooh, e_o2 = row
        g_ = e_
        g_oh = e_oh + dgoh
        g_o = e_o + dgo
        g_ooh = e_ooh + dgooh 
        g_o2 = e_o2 + dgo2

        g1 = g_oh + gh2/2 - g_ - gh2o
        g2 = g_o + gh2 - g_ - gh2o
        g3 = g_ooh + gh2*3/2 - g_ - gh2o*2

        dg0 = g_o2 - g_ - go2
        dg1 = g_ooh - g_ - go2 - gh2/2 + 1.23
        dg2 = g_o + gh2o - g_ooh - gh2/2 + 1.23
        dg3 = g_oh - g_o - gh2/2 + 1.23
        dg4 = g_ + gh2o - g_oh - gh2/2 + 1.23

        dg14 = [dg1, dg2, dg3, dg4]
        dg04 = [dg0, dg1-dg0, dg2, dg3, dg4]
        overpotential14 = max(dg14)
        overpotential04 = max(dg04)
        
        # CSV에 데이터 작성 (소숫점 모두 포함)
        writer.writerow({
            'catal': catal,
            'goh': g1,
            'go': g2,
            'gooh': g3,
            'dg0': dg0,
            'dg1': dg1,
            'dg2': dg2,
            'dg3': dg3,
            'dg4': dg4,
            'op14': overpotential14,
            'op04': overpotential04
        })

# TSV 파일로 저장 (소숫점 두자리까지만)
with open('orr_overpotential_results.tsv', 'w', newline='') as tsvfile:
    fieldnames = ['catal', 'goh', 'go', 'gooh', 'dg0', 'dg1', 'dg2', 'dg3', 'dg4', 'op14', 'op04']
    writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
    
    # 헤더 작성
    writer.writeheader()
    
    for row in data:
        catal, e_, e_oh, e_o, e_ooh, e_o2 = row
        g_ = e_
        g_oh = e_oh + dgoh
        g_o = e_o + dgo
        g_ooh = e_ooh + dgooh 
        g_o2 = e_o2 + dgo2

        g1 = g_oh + gh2/2 - g_ - gh2o
        g2 = g_o + gh2 - g_ - gh2o
        g3 = g_ooh + gh2*3/2 - g_ - gh2o*2

        dg0 = g_o2 - g_ - go2
        dg1 = g_ooh - g_ - go2 - gh2/2 + 1.23
        dg2 = g_o + gh2o - g_ooh - gh2/2 + 1.23
        dg3 = g_oh - g_o - gh2/2 + 1.23
        dg4 = g_ + gh2o - g_oh - gh2/2 + 1.23

        dg14 = [dg1, dg2, dg3, dg4]
        dg04 = [dg0, dg1-dg0, dg2, dg3, dg4]
        overpotential14 = max(dg14)
        overpotential04 = max(dg04)
        
        # TSV에 데이터 작성 (소숫점 두자리까지만)
        writer.writerow({
            'catal': catal,
            'goh': round(g1, 2),
            'go': round(g2, 2),
            'gooh': round(g3, 2),
            'dg0': round(dg0, 2),
            'dg1': round(dg1, 2),
            'dg2': round(dg2, 2),
            'dg3': round(dg3, 2),
            'dg4': round(dg4, 2),
            'op14': round(overpotential14, 2),
            'op04': round(overpotential04, 2)
        })

print("결과가 'orr_overpotential_results.csv' (소숫점 모두 포함)와 'orr_overpotential_results.tsv' (소숫점 두자리) 파일에 저장되었습니다.")

# goh와 dg0 간의 선형 관계 플롯
def plot_scaling_relationship():
    # 데이터 수집
    catalysts = []
    goh_values = []
    dg0_values = []
    
    for row in data:
        catal, e_, e_oh, e_o, e_ooh, e_o2 = row
        g_ = e_
        g_oh = e_oh + dgoh
        g_o = e_o + dgo
        g_ooh = e_ooh + dgooh 
        g_o2 = e_o2 + dgo2

        g1 = g_oh + gh2/2 - g_ - gh2o
        g2 = g_o + gh2 - g_ - gh2o
        g3 = g_ooh + gh2*3/2 - g_ - gh2o*2

        dg0 = g_o2 - g_ - go2
        
        catalysts.append(catal)
        goh_values.append(g1)
        dg0_values.append(dg0)
    
    # 플롯 생성 (scaling.py와 동일한 스타일)
    plt.figure(figsize=(4, 3), dpi=300)
    plt.scatter(goh_values, dg0_values, s=20, zorder=4, color='blue')
    
    # 촉매 이름 표시
    for i, (catal, goh, dg0) in enumerate(zip(catalysts, goh_values, dg0_values)):
        plt.annotate(f'{catal}', (goh, dg0), textcoords="offset points", 
                     xytext=(0, 5), ha='center', color='black', fontsize=8, zorder=5)
    
    # 선형 회귀
    coeffs = np.polyfit(goh_values, dg0_values, 1)
    line = np.poly1d(coeffs)
    xx = np.linspace(min(goh_values), max(goh_values), 100)
    plt.plot(xx, line(xx), label=fr'$\Delta$G$_{{\sf O_2}}$ (trend)', linestyle='-', color='black', zorder=3)
    
    # R² 계산
    y_pred = line(goh_values)
    ss_total = np.sum((np.array(dg0_values) - np.mean(dg0_values))**2)
    ss_residual = np.sum((np.array(dg0_values) - y_pred)**2)
    r_squared = 1 - (ss_residual / ss_total)
    
    # 방정식과 R² 표시
    equation = f'y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}'
    equation_with_r2 = f'{equation} (R² = {r_squared:.2f})'
    plt.text(0.34, 0.11, equation_with_r2, transform=plt.gca().transAxes, fontsize=10, color='black')
    
    plt.xlabel(fr'$\Delta$G$_{{\sf OH}}$ (eV)', fontsize='large')
    plt.ylabel(fr'$\Delta$G$_{{\sf O_2}}$ (eV)', fontsize='large')
    plt.xlim(0.1, 1.3)
    plt.ylim(-0.9, 0.0)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tight_layout()
    plt.savefig('orr_scaling_relationship.png')
    print("Scaling relationship plot saved as 'orr_scaling_relationship.png'")
    plt.close()

# 플롯 함수 실행
plot_scaling_relationship()