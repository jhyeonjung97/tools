#!/usr/bin/env python3
"""ICOHPLIST: COHP#별(sd×p) 합을 금속 원자 단위로 다시 합산.

icohp3.py와 동일한 결합별(sd 금속 × p 산소) ICOHP를 쓰되,
산소 번호가 달라도 같은 금속(Ru1 등)이면 모두 더한다.
(예: COHP#1–5의 Ru1–O19/O20/… 합 → Ru1 한 줄)
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import icohp3 as ich


def coarse_metal_label(am: str, an: str) -> str | None:
    """coarse 쌍에서 금속 쪽 라벨 (Ru1, Hf16 …)."""
    for lab in (am, an):
        m = ich.ATOM_MU_RE.match(lab)
        if m and m.group(1) != "O":
            return lab
    return None


def sort_key_by_idx(lab: str) -> tuple:
    """출력 순서: idx 오름차순, 같으면 원소 기호."""
    ele, idx_s = ich.coarse_label_ele_idx(lab)
    if idx_s == "?":
        return (1 << 30, ele)
    try:
        return (int(idx_s), ele)
    except ValueError:
        return (1 << 30, ele)


def main() -> None:
    p = argparse.ArgumentParser(
        description="금속 원자별 sum: 각 COHP#의 (sd×p) 합을 금속으로 모음"
    )
    p.add_argument(
        "icohplist",
        type=Path,
        nargs="?",
        default="ICOHPLIST.lobster",
        help="ICOHPLIST.lobster 경로",
    )
    p.add_argument(
        "--spin",
        type=int,
        choices=(1, 2),
        default=None,
        help="해당 스핀만 (미지정 시 스핀 1+2)",
    )
    p.add_argument(
        "--idx1",
        type=int,
        default=None,
        metavar="N",
        help="icohp3와 같이 coarse atomMU(am0) 원자 번호가 N인 COHP만 집계",
    )
    args = p.parse_args()

    if not args.icohplist.is_file():
        raise SystemExit(f"파일 없음: {args.icohplist}")

    spin_note = f" (spin {args.spin}만)" if args.spin else " (spin 1+2)"
    cohp_sums, meta = ich.sum_by_cohp_sd_mu_p_nu(args.icohplist, args.spin)

    by_metal: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    dist_sum: dict[str, float] = defaultdict(float)
    for cid, val in cohp_sums.items():
        if cid not in meta:
            continue
        am0, an0, dist = meta[cid]
        if args.idx1 is not None:
            i1 = ich.idx1_from_atom_mu_label(am0)
            if i1 != args.idx1:
                continue
        mlab = coarse_metal_label(am0, an0)
        if not mlab:
            continue
        by_metal[mlab] += val
        counts[mlab] += 1
        dist_sum[mlab] += dist

    # print(f"# {args.icohplist}{spin_note}")
    # print(
    #     "# 금속별 합: icohp3와 동일한 COHP#별 (금속 s·d × 산소 p)를 "
    #     "같은 금속 원자로 합산 (산소 번호 무관). d_mean = 해당 결합들의 distance 평균. "
    #     "행 순서 = idx 오름차순"
    # )
    print(
        f"# {'ele':>4} {'idx':>5}  n_bonds  {'d_mean':>10}  sum_ICOHP(sd×p)"
    )
    for lab in sorted(by_metal, key=sort_key_by_idx):
        ele, idx = ich.coarse_label_ele_idx(lab)
        n = counts[lab]
        d_mean = dist_sum[lab] / n if n else 0.0
        print(
            f"  {ele:>4} {idx:>5}  {n:7d}  {d_mean:10.5f}  {by_metal[lab]:.5f}"
        )


if __name__ == "__main__":
    main()
