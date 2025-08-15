#!/bin/env python3
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D

data = [
    # [catal, e_, e_oh, e_o, e_ooh, e_o2]
    ['Mn', -281.649211, -292.610775, -287.435762, -296.919792, -292.465621],
    ['Fe', -280.176972, -290.947634, -285.944777, -295.188869, -290.615657],
    ['Co', -279.547869, -289.617564, -284.794233, -294.191459, -289.898750],
    ['Ni', -277.208525, -287.284249, -281.314602, -291.717757, -287.374188],
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

# print(go2)

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

# 결과를 저장할 리스트
results = []

print('catal\tgoh\tgo\tgooh\tdg0\tdg1\tdg2\tdg3\tdg4\top14\top04')
for row in data:
    catal, e_, e_oh, e_o, e_ooh, e_o2 = row
    g_ = e_
    g_oh = e_oh + dgoh
    g_o = e_o + dgo
    g_ooh = e_ooh + dgooh 
    g_o2 = e_o2 + dgo2

    g1 = g_oh + gh2/2 - g_ - gh2o # goh
    g2 = g_o + gh2 - g_ - gh2o # go
    g3 = g_ooh + gh2*3/2 - g_ - gh2o*2 # gooh

    dg0 = g_o2 - g_ - go2 # go2
    dg1 = g_ooh - g_ - go2 - gh2/2 + 1.23
    dg2 = g_o + gh2o - g_ooh - gh2/2 + 1.23
    dg3 = g_oh - g_o - gh2/2 + 1.23
    dg4 = g_ + gh2o - g_oh - gh2/2 + 1.23

    dg14 = [dg1, dg2, dg3, dg4]
    dg04 = [dg0, dg1-dg0, dg2, dg3, dg4]
    overpotential14 = max(dg14)
    overpotential04 = max(dg04)
    
    # 결과를 리스트에 저장
    result_row = [
        catal, 
        round(g1, 2), 
        round(g2, 2), 
        round(g3, 2), 
        round(dg0, 2), 
        round(dg1, 2), 
        round(dg2, 2), 
        round(dg3, 2), 
        round(dg4, 2), 
        round(overpotential14, 2), 
        round(overpotential04, 2)
    ]
    results.append(result_row)
    
    print(f'{catal}\t{round(g1, 2)}\t{round(g2, 2)}\t{round(g3, 2)}\t{round(dg0, 2)}\t{round(dg1, 2)}\t{round(dg2, 2)}\t{round(dg3, 2)}\t{round(dg4, 2)}\t{round(overpotential14, 2)}\t{round(overpotential04, 2)}')

# CSV 파일로 저장
headers = ['catal', 'goh', 'go', 'gooh', 'dg0', 'dg1', 'dg2', 'dg3', 'dg4', 'op14', 'op04']

with open('orr_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(headers)
    writer.writerows(results)

# TSV 파일로 저장 (소숫점 두자리)
with open('orr_results.tsv', 'w', newline='', encoding='utf-8') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(headers)
    writer.writerows(results)

print(f"\n결과가 'orr_results.csv'와 'orr_results.tsv' 파일로 저장되었습니다.")

# goh와 dg0 간의 선형관계 그래프 그리기
def plot_scaling_relationship():
    # 데이터 추출
    goh_values = [result[1] for result in results]  # goh 값
    dg0_values = [result[4] for result in results]  # dg0 값
    catalysts = [result[0] for result in results]   # 촉매 이름
    
    # 선형 회귀
    coeffs = np.polyfit(goh_values, dg0_values, 1)
    line = np.poly1d(coeffs)
    
    # 그래프 그리기
    plt.figure(figsize=(4, 3), dpi=300)
    
    # 산점도
    plt.scatter(goh_values, dg0_values, c='blue', s=20, zorder=4)
    
    # 촉매 이름 라벨
    for xi, yi, metal in zip(goh_values, dg0_values, catalysts):
        plt.annotate(f'{metal}', (xi, yi), textcoords="offset points", 
                     xytext=(0, 5), ha='center', color='black', fontsize=8, zorder=5)
    
    # 추세선
    xx = np.linspace(min(goh_values), max(goh_values), 100)
    plt.plot(xx, line(xx), label=fr'$\Delta$G$_{{0}}$ (trend)', linestyle='-', color='black', zorder=3)
    
    # 방정식과 R² 계산
    equation = f'y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}'
    y_pred = line(goh_values)
    ss_total = np.sum((dg0_values - np.mean(dg0_values))**2)
    ss_residual = np.sum((dg0_values - y_pred)**2)
    r_squared = 1 - (ss_residual / ss_total)
    equation_with_r2 = f'{equation} (R² = {r_squared:.2f})'
    
    plt.text(0.30, 0.1, equation_with_r2, transform=plt.gca().transAxes, fontsize=10, color='black')
    
    plt.xlim(0.1, 1.3)
    plt.ylim(-0.9, 0.0)
    # 축 라벨
    plt.xlabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize='large')
    plt.ylabel(r'$\Delta$G$_{\sf OO}$ (eV)', fontsize='large')
    
    # 축 포맷
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    plt.tight_layout()
    plt.savefig('orr_goh_dg0_scaling.png')
    print("Figure saved as orr_goh_dg0_scaling.png")
    plt.close()

# g3와 dg0 간의 선형관계 그래프 그리기
def plot_g3_dg0_scaling():
    # 데이터 추출
    g3_values = [result[3] for result in results]  # g3 (gooh) 값
    dg0_values = [result[4] for result in results]  # dg0 값
    catalysts = [result[0] for result in results]   # 촉매 이름
    
    # 선형 회귀
    coeffs = np.polyfit(g3_values, dg0_values, 1)
    line = np.poly1d(coeffs)
    
    # 그래프 그리기
    plt.figure(figsize=(4, 3), dpi=300)
    
    # 산점도
    plt.scatter(g3_values, dg0_values, c='blue', s=20, zorder=4)
    
    # 촉매 이름 라벨
    for xi, yi, metal in zip(g3_values, dg0_values, catalysts):
        plt.annotate(f'{metal}', (xi, yi), textcoords="offset points", 
                     xytext=(0, 5), ha='center', color='black', fontsize=8, zorder=5)
    
    # 추세선
    xx = np.linspace(min(g3_values), max(g3_values), 100)
    plt.plot(xx, line(xx), label=fr'$\Delta$G$_{{0}}$ (trend)', linestyle='-', color='black', zorder=3)
    
    # 방정식과 R² 계산
    equation = f'y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}'
    y_pred = line(g3_values)
    ss_total = np.sum((dg0_values - np.mean(dg0_values))**2)
    ss_residual = np.sum((dg0_values - y_pred)**2)
    r_squared = 1 - (ss_residual / ss_total)
    equation_with_r2 = f'{equation} (R² = {r_squared:.2f})'
    
    plt.text(0.30, 0.1, equation_with_r2, transform=plt.gca().transAxes, fontsize=10, color='black')
    plt.xlim(3.3, 4.3)
    plt.ylim(-0.9, 0.0)
    # 축 라벨
    plt.xlabel(r'$\Delta$G$_{\sf OOH}$ (eV)', fontsize='large')
    plt.ylabel(r'$\Delta$G$_{\sf OO}$ (eV)', fontsize='large')
    
    # 축 포맷
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    plt.tight_layout()
    plt.savefig('orr_g3_dg0_scaling.png')
    print("Figure saved as orr_g3_dg0_scaling.png")
    plt.close()

# 그래프 그리기 실행
plot_scaling_relationship()
plot_g3_dg0_scaling()