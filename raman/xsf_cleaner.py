import math
import argparse
from pathlib import Path


def clean_xsf(input_path: str, cutoff: float, output_prefix: str = "clean_") -> Path:
    """
    XSF(특히 PRIMCOORD 섹션)를 읽어서, 각 원자 라인의
    변위벡터 (dx, dy, dz)의 modulus = sqrt(dx^2+dy^2+dz^2)가 cutoff 이하이면
    (dx, dy, dz) = (0, 0, 0)으로 바꾼 뒤 *_clean.xsf로 저장합니다.
    """
    input_path = str(input_path)
    path = Path(input_path)
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)

    primcoord_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "PRIMCOORD":
            primcoord_idx = i
            break
    if primcoord_idx is None:
        raise ValueError("PRIMCOORD 블록을 찾지 못했습니다.")
    if primcoord_idx + 2 >= len(lines):
        raise ValueError("PRIMCOORD 헤더/원자 라인이 충분하지 않습니다.")

    header_tokens = lines[primcoord_idx + 1].split()
    if len(header_tokens) < 2:
        raise ValueError(
            "PRIMCOORD 헤더 라인이 'N M' 형태가 아닙니다: " + lines[primcoord_idx + 1].strip()
        )
    n_atoms = int(header_tokens[0])

    atom_start = primcoord_idx + 2
    atom_end = atom_start + n_atoms
    if atom_end > len(lines):
        raise ValueError(
            f"예상 원자 수({n_atoms})만큼 라인이 없습니다. 현재 lines={len(lines)}"
        )

    kept = 0
    new_lines = list(lines)

    for j in range(n_atoms):
        line_idx = atom_start + j
        toks = new_lines[line_idx].split()
        if len(toks) < 7:
            raise ValueError(
                f"원자 라인 토큰 수가 7 미만입니다: line={new_lines[line_idx].rstrip()}"
            )

        elem = toks[0]
        # x y z dx dy dz
        x, y, z, dx, dy, dz = map(float, toks[1:7])
        mod = math.sqrt(dx * dx + dy * dy + dz * dz)

        if mod <= cutoff:
            dx2 = dy2 = dz2 = 0.0
        else:
            dx2, dy2, dz2 = dx, dy, dz
            kept += 1

        # float를 고정밀로 출력
        new_lines[line_idx] = (
            f"{elem:2s} "
            f"{x: .16f} {y: .16f} {z: .16f} "
            f"{dx2: .16f} {dy2: .16f} {dz2: .16f}\n"
        )

    if output_prefix:
        out_path = path.with_name(output_prefix + path.name)
    else:
        out_path = path

    out_path.write_text("".join(new_lines), encoding="utf-8")

    print(f"입력: {path}")
    print(f"출력: {out_path}")
    print(f"cutoff: {cutoff}")
    print(f"모듈러스 > cutoff 인 원자(벡터 유지) 수: {kept} / {n_atoms}")

    return out_path


def main():
    parser = argparse.ArgumentParser(description="modulus cutoff 기준으로 XSF 변위벡터 클린")
    parser.add_argument("--input", default="mode_0008.xsf", help="입력 XSF 파일 경로")
    parser.add_argument("--cutoff", type=float, default=0.05, help="modulus가 cutoff 이하이면 (0,0,0)으로")
    parser.add_argument("--all", action="store_true", help="dir/pattern 아래 모든 mode_*.xsf 처리")
    parser.add_argument("--dir", default=".", help="--all 사용 시 파일을 찾는 디렉토리")
    parser.add_argument("--pattern", default="mode_*.xsf", help="--all 사용 시 글롭 패턴")
    parser.add_argument("--prefix", default="clean_", help="출력 파일 prefix (예: clean_)")
    args = parser.parse_args()

    if args.all:
        base_dir = Path(args.dir)
        paths = sorted(base_dir.glob(args.pattern))
        if not paths:
            raise FileNotFoundError(f"--all 활성화 후 패턴 '{args.pattern}'에 해당하는 파일이 없습니다: {base_dir}")
        for p in paths:
            clean_xsf(str(p), args.cutoff, output_prefix=args.prefix)
    else:
        clean_xsf(args.input, args.cutoff, output_prefix=args.prefix)


if __name__ == "__main__":
    main()

