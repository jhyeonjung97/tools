from ase.io import read, write
from ase.constraints import FixAtoms

# POSCAR 파일 읽기
atoms = read('POSCAR')

# N, Co 원소에 해당하는 원자들의 z 방향 고정
mask = [atom.symbol in ['N', 'Co'] for atom in atoms]  # N과 Co 원소에 대해 True 설정

# ASE의 constraints 적용
fix_atoms = FixAtoms(mask=mask)  # 선택된 원자들만 고정
atoms.set_constraint(fix_atoms)

# 수정된 구조를 a.vasp로 저장 (VASP 5 포맷)
write('a.vasp', atoms, format='vasp', vasp5=True)

print("a.vasp 파일이 성공적으로 생성되었습니다.")