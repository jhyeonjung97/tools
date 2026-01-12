#!/usr/bin/env python3
"""
현재 폴더 내의 모든 JSON 파일을 CIF 형식으로 변환하는 스크립트
"""

from ase.io import read, write
import os
import glob
import argparse

def convert_json_to_cif(json_file, output_dir=None, overwrite=False):
    """JSON 파일을 CIF로 변환하는 함수
    
    Parameters:
    -----------
    json_file : str
        입력 JSON 파일 경로
    output_dir : str, optional
        출력 디렉토리 (None이면 입력 파일과 같은 디렉토리)
    overwrite : bool, optional
        기존 파일을 덮어쓸지 여부
    """
    try:
        # JSON 파일 읽기
        atoms = read(json_file)
        
        # 출력 파일명 생성
        base_name = os.path.splitext(os.path.basename(json_file))[0]
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            cif_file = os.path.join(output_dir, f"{base_name}.cif")
        else:
            cif_file = json_file.replace('.json', '.cif')
        
        # 파일이 이미 존재하는 경우 확인
        if os.path.exists(cif_file) and not overwrite:
            print(f"Skip: {cif_file} already exists (use --overwrite to overwrite)")
            return False
        
        # CIF 파일로 저장
        write(cif_file, atoms, format='cif')
        print(f"Converted: {json_file} -> {cif_file}")
        return True
        
    except Exception as e:
        print(f"Error converting {json_file}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='현재 폴더 내의 모든 JSON 파일을 CIF 형식으로 변환',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  python json2cif.py                    # 현재 폴더의 모든 JSON 파일 변환
  python json2cif.py --output ./cif     # 출력을 ./cif 폴더에 저장
  python json2cif.py --overwrite        # 기존 CIF 파일 덮어쓰기
  python json2cif.py --pattern "*.json" # 특정 패턴의 파일만 변환
        """
    )
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='출력 디렉토리 (기본값: 입력 파일과 같은 디렉토리)')
    parser.add_argument('--overwrite', action='store_true',
                        help='기존 CIF 파일을 덮어쓰기')
    parser.add_argument('--pattern', type=str, default='*.json',
                        help='변환할 파일 패턴 (기본값: *.json)')
    parser.add_argument('--directory', '-d', type=str, default='.',
                        help='작업할 디렉토리 (기본값: 현재 디렉토리)')
    
    args = parser.parse_args()
    
    # JSON 파일 찾기
    pattern = os.path.join(args.directory, args.pattern)
    json_files = glob.glob(pattern)
    
    if not json_files:
        print(f"No JSON files found matching pattern: {pattern}")
        return
    
    print(f"Found {len(json_files)} JSON file(s)")
    print("-" * 60)
    
    # 각 파일 변환
    success_count = 0
    for json_file in sorted(json_files):
        if convert_json_to_cif(json_file, args.output, args.overwrite):
            success_count += 1
    
    print("-" * 60)
    print(f"Conversion complete: {success_count}/{len(json_files)} files converted successfully")

if __name__ == '__main__':
    main()




