"""
Pourbaix 형태의 pH–전위도: x = pH, y = 전위 (V vs. SHE).
물의 안정 영역(HER/OER 평형선)을 예시로 표시합니다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

# 298.15 K에서 RT/F * ln(10) [V] (Nernst 기울기)
CONST = 8.617e-5 * 298.15 * np.log(10)


def water_stability_vs_she(pH: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """2H⁺ + 2e⁻ = H₂, O₂ + 4H⁺ + 4e⁻ = 2H₂O 평형선 (대략 1.23 V 기준)."""
    e_her = 0.0 - CONST * pH
    e_oer = 1.23 - CONST * pH
    return e_her, e_oer


def main() -> None:
    pH = np.linspace(0, 14, 400)
    e_her, e_oer = water_stability_vs_she(pH)
    u_min, u_max = -1.0, 2.0

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.set_xlim(0, 14)
    ax.set_ylim(u_min, u_max)
    ax.set_xlabel("pH")
    ax.set_ylabel("Potential (V vs. SHE)")

    ax.fill_between(pH, e_her, e_oer, color="silver", alpha=0.5)
    ax.plot(pH, e_her, color="blue", ls="-", lw=1.0, label="HER")
    ax.plot(pH, e_oer, color="red", ls="-", lw=1.0, label="OER")

    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)

    plt.tight_layout()
    plt.savefig("mat303_pourbaix.png", dpi=220, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
