#!/usr/bin/env python3
"""
RuO2 시스템의 Bader charge 분석 스크립트
0_, 1_, 2_, 3_, 4_ 폴더에서 atoms_bader_charge.json 파일을 읽어서
각 원소별 평균 charge를 계산하고 표로 출력합니다.
"""

import os
import json
import numpy as np
import pandas as pd
from ase.io import read

def get_element_symbol(atomic_number):
    """원자 번호를 원소 기호로 변환"""
    element_map = {44: 'Ru', 75: 'Re', 8: 'O'}
    return element_map.get(atomic_number, f'Unknown_{atomic_number}')

def analyze_bader_charges(base_path):
    """Bader charge 데이터를 분석하여 원소별 평균 charge를 계산"""
    
    # 결과를 저장할 딕셔너리
    results = {}
    
    # 각 폴더 (0_, 1_, 2_, 3_, 4_) 처리
    for folder in ['0_', '1_', '2_', '3_', '4_']:
        json_path = os.path.join(base_path, folder, 'atoms_bader_charge.json')
        
        if not os.path.exists(json_path):
            print(f"경고: {json_path} 파일을 찾을 수 없습니다.")
            continue
            
        try:
            # ASE로 JSON 파일 읽기
            atoms = read(json_path)
            
            # initial_charges와 numbers 추출
            # ASE가 __ndarray__ 구조를 자동으로 파싱해줍니다
            initial_charges = atoms.get_initial_charges()
            numbers = atoms.get_atomic_numbers()
            
            if len(initial_charges) != len(numbers):
                print(f"경고: {folder}에서 charge({len(initial_charges)})와 atomic number({len(numbers)}) 개수가 일치하지 않습니다.")
                continue
            
            # 원소별로 charge 분류
            element_charges = {}
            for charge, atomic_num in zip(initial_charges, numbers):
                element = get_element_symbol(atomic_num)
                if element not in element_charges:
                    element_charges[element] = []
                element_charges[element].append(charge)
            
            # 결과 저장
            results[folder] = element_charges
            
            
        except Exception as e:
            print(f"오류: {folder} 폴더 처리 중 오류 발생: {e}")
            continue
    
    return results

def create_summary_table(results):
    """결과를 표 형태로 정리"""
    
    # 모든 원소 수집
    all_elements = set()
    for folder_data in results.values():
        all_elements.update(folder_data.keys())
    
    all_elements = sorted(all_elements)
    
    # 표 데이터 생성
    table_data = []
    for folder in ['0_', '1_', '2_', '3_', '4_']:
        if folder in results:
            row = {'Folder': folder}
            for element in all_elements:
                if element in results[folder]:
                    charges = results[folder][element]
                    row[f'{element}_avg'] = f"{np.mean(charges):.3f}"
                    row[f'{element}_count'] = len(charges)
                else:
                    row[f'{element}_avg'] = "N/A"
                    row[f'{element}_count'] = 0
            table_data.append(row)
    
    # DataFrame 생성
    df = pd.DataFrame(table_data)
    
    return df, all_elements

def main():
    # 현재 작업 디렉토리 사용
    base_path = os.getcwd()
    
    # Bader charge 분석
    results = analyze_bader_charges(base_path)
    
    if not results:
        return
    
    # 요약 표 생성
    df, elements = create_summary_table(results)
    
    # 표 출력 (평균 Charge만)
    avg_cols = ['Folder'] + [f'{elem}_avg' for elem in elements]
    print(df[avg_cols].to_string(index=False))

if __name__ == "__main__":
    main()
