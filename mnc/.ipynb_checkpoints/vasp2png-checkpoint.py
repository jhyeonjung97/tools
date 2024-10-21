import os
from ase.io import read, write

# 현재 디렉토리의 모든 .vasp 파일 읽기
vasp_files = [f for f in os.listdir() if f.endswith('.vasp')]

# 각 VASP 파일을 PNG로 저장
for vasp_file in vasp_files:
    # VASP 파일 읽기
    atoms = read(vasp_file)
    
    # PNG 파일명 설정
    png_file = vasp_file.replace('.vasp', '.png')
    
    # ASE의 write 함수를 사용해 x축에서 본 구조를 PNG로 저장
    write(png_file, atoms, rotation='0x, 0y, 0z')

    print(f"{png_file}로 저장되었습니다.")
