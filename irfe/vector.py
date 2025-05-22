#!/usr/bin/env python
import os
import glob
import numpy as np
import shutil
from ase.io import read

def compare_atomic_positions(base_path):
    """
    지정된 경로 패턴에 맞는 폴더에서 CONTCAR와 vib/POSCAR 파일의 원자 위치를 비교합니다.
    위치가 다른 경우 해당 경로를 출력합니다.
    
    Args:
        base_path: 검색할 기본 경로 패턴
    """
    # 지정된 패턴에 맞는 모든 폴더 찾기
    all_dirs = glob.glob(base_path)
    
    for dir_path in all_dirs:
        done_path = os.path.join(dir_path, 'DONE')
        contcar_path = os.path.join(dir_path, 'CONTCAR')
        poscar_path = os.path.join(dir_path, 'vib', 'POSCAR')
        
        # DONE 파일이 있는지 확인
        if not os.path.exists(done_path):
            continue
        
        # 파일 존재 여부 확인
        contcar_exists = os.path.exists(contcar_path)
        poscar_exists = os.path.exists(poscar_path)
        
        if not contcar_exists and not poscar_exists:
            print(f"두 파일 모두 없음: {dir_path}")
            continue
        elif not contcar_exists:
            print(f"CONTCAR 파일 없음: {dir_path}")
            continue
        elif not poscar_exists:
            print(f"vib/POSCAR 파일 없음: {dir_path}")
            continue
        
        try:
            # ASE를 사용하여 파일 읽기
            contcar = read(contcar_path, format='vasp')
            poscar = read(poscar_path, format='vasp')
            
            # 원자 위치 비교
            contcar_positions = contcar.get_positions()
            poscar_positions = poscar.get_positions()
            
            # 원자 수가 같은지 확인
            if len(contcar_positions) != len(poscar_positions):
                print(f"원자 수가 다릅니다: {dir_path}")
                continue
            
            # 위치 차이 확인 (작은 수치 오차 허용)
            if not np.allclose(contcar_positions, poscar_positions, atol=1e-5):
                print(f"원자 위치가 다릅니다: {dir_path}")
                vib_folder = os.path.join(dir_path, 'vib')
                if os.path.exists(vib_folder):
                    try:
                        shutil.rmtree(vib_folder)
                        print(f"vib 폴더 삭제 완료: {vib_folder}")
                    except Exception as e:
                        print(f"vib 폴더 삭제 중 오류 발생: {vib_folder} - {str(e)}")
        
        except Exception as e:
            print(f"파일 처리 중 오류 발생: {dir_path} - {str(e)}")

if __name__ == "__main__":
    # 사용자가 지정한 경로 패턴
    path_pattern = "/home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*"
    compare_atomic_positions(path_pattern)
    print("비교 완료!")
