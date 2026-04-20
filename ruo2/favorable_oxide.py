import os
import sys
import mendeleev as md
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
from contextlib import nullcontext
import warnings

# MPRester import - 최신 버전과 구버전 모두 지원
try:
    from mp_api.client import MPRester
except ImportError:
    try:
        from pymatgen.ext.matproj import MPRester
    except ImportError:
        raise ImportError("MPRester를 찾을 수 없습니다. mp-api 또는 pymatgen을 설치하세요.")

warnings.filterwarnings('ignore')


def get_favorable_oxide(element, api_key=None):
    """
    Materials Project에서 해당 element의 oxide 중 가장 낮은 formation energy를 가진 구조를 찾습니다.
    
    Args:
        element: 원소 기호 (예: 'Ti', 'V')
        api_key: Materials Project API 키 (없으면 환경변수에서 가져옴)
        
    Returns:
        dict: 가장 낮은 formation energy를 가진 oxide 정보
            {'formula': 'TiO2', 'formation_energy': -2.5, 'entry': entry_object, 'material_id': 'mp-1234'}
    """
    if api_key is None:
        api_key = os.getenv('MAPI_KEY')
        if not api_key:
            print(f"경고: {element}에 대해 MAPI_KEY 환경변수가 설정되지 않았습니다.")
            return None
    
    print(f"{element}-O 시스템의 데이터를 가져오는 중...")
    
    # Materials Project에서 entries 가져오기
    entries = None
    use_old_api = False
    
    try:
        with MPRester(api_key) as mpr:
            entries = mpr.get_entries_in_chemsys([element, 'O'])
            print(f"  총 {len(entries)}개의 entries를 가져왔습니다.")
    except AttributeError as e:
        if "model_fields" in str(e):
            print("  mp_api 라이브러리 버전 호환성 문제. 구버전 API를 시도합니다...")
            use_old_api = True
        else:
            raise
    except Exception as e:
        print(f"  에러 발생: {type(e).__name__}: {e}")
        print("  구버전 API를 시도합니다...")
        use_old_api = True
    
    # 구버전 API 사용
    if use_old_api or entries is None:
        try:
            from pymatgen.ext.matproj import MPRester as OldMPRester
            with OldMPRester(api_key) as old_mpr:
                entries = old_mpr.get_entries_in_chemsys([element, 'O'])
                print(f"  총 {len(entries)}개의 entries를 가져왔습니다 (구버전 API 사용).")
        except Exception as e2:
            print(f"  구버전 API도 실패: {e2}")
            return None
    
    if entries is None or len(entries) == 0:
        print(f"  {element}-O 시스템의 데이터를 찾을 수 없습니다.")
        return None
    
    # Phase diagram 생성
    try:
        pd = PhaseDiagram(entries)
    except Exception as e:
        print(f"  Phase diagram 생성 실패: {e}")
        return None
    
    # Oxide만 필터링 (O를 포함하고 해당 element를 포함하는 entries)
    oxide_entries = []
    for entry in entries:
        comp = entry.composition
        if 'O' in comp and element in comp:
            # Pure element는 제외
            if len(comp.elements) > 1:
                oxide_entries.append(entry)
    
    if len(oxide_entries) == 0:
        print(f"  {element}의 oxide를 찾을 수 없습니다.")
        return None
    
    # 각 oxide의 formation energy 계산
    oxide_data = []
    for entry in oxide_entries:
        try:
            # Phase diagram을 기준으로 formation energy 계산
            comp = entry.composition
            formula = comp.reduced_formula
            
            # Formation energy 계산 (여러 방법 시도)
            try:
                # 방법 1: get_form_energy를 사용하고 원자 수로 나누기
                formation_energy = pd.get_form_energy(entry) / comp.num_atoms
            except Exception as e1:
                try:
                    # 방법 2: get_form_energy_per_atom 시도 (Composition 객체로)
                    formation_energy = pd.get_form_energy_per_atom(comp)
                except Exception as e2:
                    # 방법 3: entry의 energy_per_atom 사용 (formation energy가 아님, 경고)
                    print(f"  경고: {formula}의 formation energy 계산 실패, energy_per_atom 사용")
                    formation_energy = entry.energy_per_atom
            
            # Material ID 추출
            material_id = None
            # 여러 방법으로 material_id 추출 시도
            if hasattr(entry, 'material_id'):
                material_id = entry.material_id
            elif hasattr(entry, 'entry_id'):
                entry_id = str(entry.entry_id)
                if 'mp-' in entry_id:
                    # mp-1234 형식 추출
                    parts = entry_id.split('mp-')
                    if len(parts) > 1:
                        mp_num = parts[1].split('-')[0].split('_')[0]
                        material_id = f'mp-{mp_num}'
            elif hasattr(entry, 'structure'):
                # structure에서 material_id 가져오기 시도
                struct = entry.structure
                if hasattr(struct, 'material_id'):
                    material_id = struct.material_id
            
            # 여전히 없으면 entry의 name이나 다른 속성에서 추출 시도
            if material_id is None:
                entry_str = str(entry)
                if 'mp-' in entry_str:
                    import re
                    match = re.search(r'mp-\d+', entry_str)
                    if match:
                        material_id = match.group()
            
            oxide_data.append({
                'formula': formula,
                'formation_energy': formation_energy,
                'entry': entry,
                'material_id': material_id,
                'energy_per_atom': entry.energy_per_atom
            })
        except Exception as e:
            # 에러 메시지는 조용하게 처리 (너무 많은 출력 방지)
            # print(f"  {entry.composition.reduced_formula}의 formation energy 계산 실패: {e}")
            continue
    
    if len(oxide_data) == 0:
        print(f"  {element}의 oxide formation energy를 계산할 수 없습니다.")
        return None
    
    # Formation energy가 가장 낮은 oxide 찾기
    best_oxide = min(oxide_data, key=lambda x: x['formation_energy'])
    
    print(f"  가장 낮은 formation energy: {best_oxide['formula']} "
          f"({best_oxide['formation_energy']:.4f} eV/atom)")
    
    return best_oxide


def download_structure(material_id, api_key, output_dir, formula=None):
    """
    Materials Project에서 구조를 다운로드하고 CIF 파일로 저장합니다.
    
    Args:
        material_id: Materials Project material ID (예: 'mp-1234')
        api_key: Materials Project API 키
        output_dir: 저장할 디렉토리 경로
        formula: 화학식 (파일명에 사용, None이면 material_id 사용)
        
    Returns:
        str: 저장된 파일 경로, 실패 시 None
    """
    if material_id is None:
        print("  Material ID가 없어 구조를 다운로드할 수 없습니다.")
        return None
    
    # Material ID 정규화 (mp- 접두사 확인)
    if not material_id.startswith('mp-'):
        material_id = f'mp-{material_id}'
    
    print(f"  구조 다운로드 중: {material_id}...")
    
    adaptor = AseAtomsAdaptor()
    use_old_api = False
    
    try:
        with MPRester(api_key) as mpr:
            try:
                structure = mpr.get_structure_by_material_id(material_id)
            except AttributeError:
                # 구버전 API 시도
                use_old_api = True
                raise
    
    except (AttributeError, Exception) as e:
        if not use_old_api:
            print(f"  최신 API 실패, 구버전 API 시도 중...")
            use_old_api = True
        
        if use_old_api:
            try:
                from pymatgen.ext.matproj import MPRester as OldMPRester
                with OldMPRester(api_key) as old_mpr:
                    structure = old_mpr.get_structure_by_material_id(material_id)
            except Exception as e2:
                print(f"  구조 다운로드 실패: {e2}")
                return None
    
    if structure is None:
        print(f"  {material_id}에 대한 구조를 찾을 수 없습니다.")
        return None
    
    # pymatgen Structure를 ASE Atoms로 변환
    try:
        atoms = adaptor.get_atoms(structure)
    except Exception as e:
        print(f"  구조 변환 실패: {e}")
        return None
    
    # 파일명 생성
    if formula:
        # 화학식에서 특수문자 제거 (예: TiO2 -> TiO2)
        safe_formula = formula.replace(' ', '').replace('(', '').replace(')', '')
        filename = f"{safe_formula}_{material_id}.cif"
    else:
        filename = f"{material_id}.cif"
    
    output_path = os.path.join(output_dir, filename)
    
    # CIF 파일로 저장
    try:
        write(output_path, atoms, format='cif')
        print(f"  구조가 {output_path}에 저장되었습니다.")
        return output_path
    except Exception as e:
        print(f"  파일 저장 실패: {e}")
        return None


def main():
    """메인 함수"""
    # 원소 리스트
    element_indices = list(range(21, 31)) + list(range(39, 49)) + [57] + list(range(72, 81))
    
    # API 키 확인
    api_key = os.getenv('MAPI_KEY')
    if not api_key:
        print("경고: MAPI_KEY 환경변수가 설정되지 않았습니다.")
        print("Materials Project API 키를 설정하거나 환경변수로 제공하세요.")
        response = input("계속 진행하시겠습니까? (y/n): ")
        if response.lower() != 'y':
            return
    
    # 각 element에 대해 가장 낮은 formation energy를 가진 oxide 찾기 및 다운로드
    elements_data = {}
    
    for i in element_indices:
        element = md.element(i).symbol
        print(f"\n{'='*60}")
        print(f"Processing: {element} (atomic number {i})")
        print(f"{'='*60}")
        
        oxide_data = get_favorable_oxide(element, api_key)
        elements_data[element] = oxide_data
        
        # 디렉토리 생성
        dir_name = f'{i}_{element}'
        os.makedirs(dir_name, exist_ok=True)
        
        # 결과를 파일로 저장 및 구조 다운로드
        if oxide_data:
            # 텍스트 파일로 정보 저장
            result_file = os.path.join(dir_name, 'favorable_oxide.txt')
            with open(result_file, 'w') as f:
                f.write(f"Element: {element}\n")
                f.write(f"Most Favorable Oxide: {oxide_data['formula']}\n")
                f.write(f"Formation Energy: {oxide_data['formation_energy']:.6f} eV/atom\n")
                f.write(f"Material ID: {oxide_data['material_id']}\n")
                f.write(f"Energy per Atom: {oxide_data['energy_per_atom']:.6f} eV/atom\n")
            print(f"  결과가 {result_file}에 저장되었습니다.")
            
            # 구조 다운로드
            cif_path = download_structure(
                oxide_data['material_id'], 
                api_key, 
                dir_name, 
                formula=oxide_data['formula']
            )
            if cif_path:
                print(f"  CIF 파일 저장 완료: {cif_path}")
    
    # 요약 출력
    print(f"\n{'='*60}")
    print("요약")
    print(f"{'='*60}")
    for element, oxide_data in elements_data.items():
        if oxide_data:
            print(f"{element:3s}: {oxide_data['formula']:15s} "
                  f"(E_form = {oxide_data['formation_energy']:8.4f} eV/atom, "
                  f"ID = {oxide_data['material_id']})")


if __name__ == "__main__":
    main()