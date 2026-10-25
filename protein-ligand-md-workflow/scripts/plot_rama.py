#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""画标准 Ramachandran 散点图。

为什么不用 dit xvg_show:
    gmx rama 输出的 rama.xvg 是"逐残基(或逐帧)的 (phi, psi) 点",
    dit xvg_show 只会按顺序把点连成折线, 出来是一团乱麻,
    完全失去 Ramachandran 图"看点云分布在哪些允许区"的意义。
    Ramachandran 图必须是散点。

用法:
    python3 plot_rama.py <rama.xvg> [-o analysis/rama.png] [--title "..."]

关于"允许区占比":
    背景允许区用的是常见的**矩形近似**(beta / alpha / left-handed 三块),
    不是 molprobity 那种精确等值线轮廓, 所以 in_allowed_fraction 是近似值,
    适合用来快速判断"有没有大量异常点", 不适合当作严格的统计指标写进论文。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
except Exception:
    print("错误: 需要 numpy", file=sys.stderr)
    raise SystemExit(2)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"错误: 需要 matplotlib ({exc})", file=sys.stderr)
    raise SystemExit(2)

ALLOWED_REGIONS = [
    (-180.0, -20.0, 90.0, 180.0, "beta"),
    (-160.0, -20.0, -80.0, 60.0, "alpha"),
    (-120.0, 60.0, -60.0, 0.0, "left-handed"),
]


def read_xvg(path):
    """返回 (phi, psi, labels); 只取前两列数值。"""
    phi, psi, labels = [], [], []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "@#":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                a, b = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            phi.append(a)
            psi.append(b)
            labels.append(parts[2] if len(parts) > 2 else "")
    return np.array(phi), np.array(psi), labels


def in_allowed(phi, psi):
    hit = np.zeros(len(phi), dtype=bool)
    for p0, p1, s0, s1, _ in ALLOWED_REGIONS:
        hit |= (phi >= p0) & (phi <= p1) & (psi >= s0) & (psi <= s1)
    return hit


def main(argv=None):
    ap = argparse.ArgumentParser(description="Ramachandran 散点图")
    ap.add_argument("xvg", help="gmx rama 输出的 xvg")
    ap.add_argument("-o", "--output", default=None, help="输出 png")
    ap.add_argument("--title", default="Ramachandran Plot")
    ap.add_argument("--stat", default=None, help="把统计写到这个文件")
    args = ap.parse_args(argv)

    src = Path(args.xvg)
    if not src.is_file():
        print(f"错误: 找不到 {src}", file=sys.stderr)
        return 2
    out = Path(args.output) if args.output else src.with_suffix(".png")

    phi, psi, labels = read_xvg(src)
    if len(phi) == 0:
        print(f"错误: {src} 里没有可用的数值点", file=sys.stderr)
        return 2

    frac = 100.0 * in_allowed(phi, psi).sum() / len(phi)

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    for (p0, p1, s0, s1, _), c in zip(ALLOWED_REGIONS, ["#e8f0fb", "#eef7ee", "#f7f0e8"]):
        ax.axhspan(s0, s1, color=c, zorder=0)
    ax.scatter(phi, psi, s=26, c="#1f77b4", alpha=0.65, edgecolors="none",
               zorder=3, label="residues")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.set_xlabel("Phi (deg)")
    ax.set_ylabel("Psi (deg)")
    ax.set_title(f"{args.title}\n{len(phi)} residues, {frac:.1f}% in allowed regions")
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)

    print(f"  Ramachandran 散点图 {out}  ({len(phi)} 点, 允许区占比 {frac:.1f}%)")
    if args.stat:
        Path(args.stat).write_text(
            f"points={len(phi)}\nin_allowed_fraction={frac:.2f}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())