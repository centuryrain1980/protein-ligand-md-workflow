#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""配体-蛋白界面分析与相互作用分类 (Residue Interaction Network)。

用法:
    python3 interface_analysis.py <工作目录> [选项]

选项:
    --top <file>       拓扑文件(默认依次找 md.tpr/topol.tpr/prod.tpr/complex.top)
    --traj <file>      轨迹(默认依次找 md.xtc/prod.xtc/traj.xtc)
    --out <dir>        输出目录(默认 <工作目录>/analysis)
    --lig <resname>    配体残基名(默认自动识别: 非标准残基且非水/离子)
    --prot <sel>       蛋白选择串(默认 protein)
    --cutoff 4.5       一般接触距离上限 (A)
    --hbond 2.5         氢键 H...A 距离上限 (A)
    --hbond-angle 120   氢键 D-H...A 角度下限 (度)
    --salt 4.0         盐桥距离上限 (A)
    --metal 2.8        金属配位距离上限 (A)
    --frames N         最多分析 N 帧(默认 0 = 全部)
    --stride N         抽帧步长(默认 1)

输出(写到 <out>/):
    interface_events.csv    每帧每个相互作用一行(事件级)
    interface_residues.csv  界面残基接触统计
    interface_rin.csv       残基-残基接触网络(出现频次)
    interface_summary.json  汇总统计
    interface_map.png       接触频率热图

设计原则(借鉴 PyMACS 的做法):
  * 严格区分"空间接近"与"化学上说得通的事件"——判不出来就归 Other Contact,
    不编造氢键/盐桥标签;
  * 每条事件带 ValidationStatus(confirmed / geometric / proximity-only),
    说明结论的证据强度;
  * 只用拓扑里的电荷/元素/残基名 + 几何做判断, 不做能量计算。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

try:
    import numpy as np
except Exception:
    print("错误: 需要 numpy (pip3 install numpy)", file=sys.stderr)
    raise SystemExit(2)

try:
    from scipy.spatial import cKDTree
except Exception:
    cKDTree = None

STANDARD_RES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HOH", "WAT", "SOL", "TIP3", "TIP4", "T3P", "T4P", "T5P", "SPC",
    "NA", "CL", "K", "MG", "ZN", "FE", "MN", "CU", "CA", "HSD", "HSE", "HSP",
    "ASH", "GLH", "LYN", "CYX",
}
POS_RES = {"ARG", "LYS", "HIS", "HSD", "HSE", "HSP"}
NEG_RES = {"ASP", "GLH", "ASH"}
METAL_ELEMENTS = {"NA", "MG", "ZN", "FE", "MN", "CU", "CA", "K", "CO", "NI",
                  "CD", "HG", "CU1", "CU2"}

EVENT_COLUMNS = [
    "Frame", "Time_ns", "InteractionClass", "ValidationStatus",
    "ClassificationSource", "LigandResidue", "LigandAtom", "LigandElement",
    "LigandCharge", "ProteinResidue", "ProteinAtom", "ProteinElement",
    "ProteinCharge", "Distance", "Angle", "DonorAtom", "AcceptorAtom",
]
def find_file(work, names):
    for n in names:
        p = Path(work) / n
        if p.is_file():
            return p
    return None


def norm_res(r):
    return str(r).strip().upper()


def elem_charge(atom):
    el = None
    for attr in ("element", "elements"):
        try:
            v = getattr(atom, attr)
            if v is not None:
                el = str(v).strip().upper()
                break
        except Exception:
            pass
    if not el:
        try:
            el = str(atom.name)[0].upper()
        except Exception:
            el = "X"
    ch = 0.0
    for attr in ("charge", "charges"):
        try:
            v = getattr(atom, attr)
            if v is not None:
                ch = float(v)
                break
        except Exception:
            pass
    return el, ch


def build_atom_arrays(universe, group):
    """每个原子的: 元素 / 电荷 / 名称 / 残基名 / 残基号 / 链。"""
    n = group.n_atoms
    elem = np.empty(n, dtype=object)
    chg = np.zeros(n)
    name = np.empty(n, dtype=object)
    rname = np.empty(n, dtype=object)
    rid = np.zeros(n, dtype=int)
    chain = np.empty(n, dtype=object)
    for k, a in enumerate(group.atoms):
        e, c = elem_charge(a)
        elem[k] = e
        chg[k] = c
        name[k] = a.name
        res = universe.residues[a.resindex]
        rname[k] = norm_res(res.resname)
        rid[k] = int(res.resid)
        # 某些拓扑(如 GROMACS .tpr)不带 chainID, 这时退化成用残基号做标识
        try:
            chain[k] = res.chainID
        except Exception:
            chain[k] = "A"
    return elem, chg, name, rname, rid, chain


def attached_h_map(universe):
    """原子索引 -> 与它相连的 H 的索引列表。"""
    out = {}
    try:
        bonds = universe.atoms.bonds
    except Exception:
        return out
    for b in bonds:
        i1, i2 = int(b.atoms[0].index), int(b.atoms[1].index)
        e1 = elem_charge(universe.atoms[i1])[0]
        e2 = elem_charge(universe.atoms[i2])[0]
        if e1 == "H" and e2 != "H":
            out.setdefault(i2, []).append(i1)
        elif e2 == "H" and e1 != "H":
            out.setdefault(i1, []).append(i2)
    return out


def angle_deg(a, b, c):
    """a-b-c 夹角(度), b 为顶点。"""
    v1 = a - b
    v2 = c - b
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 == 0 or n2 == 0:
        return None
    cosv = max(-1.0, min(1.0, float(np.dot(v1, v2)) / (n1 * n2)))
    return math.degrees(math.acos(cosv))


def is_donor(elem, has_h):
    return elem in ("N", "O", "S") and has_h


def is_acceptor(elem, chg, has_h):
    if elem == "O":
        return chg <= -0.4 or not has_h
    if elem == "S":
        return True
    if elem == "N":
        return (not has_h) and chg < 0.0
    return False


def salt_bridge(prn, lrn, pc, lc):
    """带强正电的残基与带强负电的残基之间。"""
    pos_side = (prn in POS_RES and pc >= 0.3) or (lrn in POS_RES and lc >= 0.3)
    neg_side = (prn in NEG_RES and pc <= -0.4) or (lrn in NEG_RES and lc <= -0.4)
    return bool(pos_side and neg_side)
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="配体-蛋白界面分析与相互作用分类 (RIN)",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("workdir", help="工作目录")
    ap.add_argument("--top", default=None)
    ap.add_argument("--traj", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--lig", default=None, help="配体残基名")
    ap.add_argument("--prot", default="protein")
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--hbond", type=float, default=2.5)
    ap.add_argument("--hbond-angle", type=float, default=120.0)
    ap.add_argument("--salt", type=float, default=4.0)
    ap.add_argument("--metal", type=float, default=2.8)
    ap.add_argument("--frames", type=int, default=0, help="0 = 全部")
    ap.add_argument("--stride", type=int, default=1)
    args = ap.parse_args(argv)

    work = Path(args.workdir).resolve()
    if not work.is_dir():
        print("错误: 工作目录不存在:", work, file=sys.stderr)
        return 2

    top = Path(args.top) if args.top else find_file(
        work, ["md.tpr", "topol.tpr", "prod.tpr", "complex.top"])
    traj = Path(args.traj) if args.traj else find_file(
        work, ["md.xtc", "prod.xtc", "traj.xtc"])
    if top is None or traj is None:
        print(f"错误: 找不到拓扑或轨迹 (top={top}, traj={traj})", file=sys.stderr)
        return 2

    try:
        import MDAnalysis as mda
    except Exception:
        print("错误: 需要 MDAnalysis (pip3 install mdanalysis)", file=sys.stderr)
        return 2

    outdir = Path(args.out) if args.out else work / "analysis"
    outdir.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(str(top), str(traj))
    prot = u.select_atoms(args.prot)
    if prot.n_atoms == 0:
        print(f"错误: 蛋白选择 '{args.prot}' 没选到原子", file=sys.stderr)
        return 2

    # ---- 配体识别 ----
    if args.lig:
        ligr = norm_res(args.lig)
    else:
        counts = Counter()
        for res in u.residues:
            rn = norm_res(res.resname)
            if rn in STANDARD_RES:
                continue
            counts[rn] += 1
        if not counts:
            print("错误: 自动识别配体失败, 请用 --lig 指定", file=sys.stderr)
            return 2
        ligr = sorted(counts.items(), key=lambda kv: kv[1])[0][0]
    lig = u.select_atoms(f"resname {ligr}")
    if lig.n_atoms == 0:
        print(f"错误: 没找到残基 {ligr}", file=sys.stderr)
        return 2

    print(f"[接口] 拓扑 {top.name} | 轨迹 {traj.name}")
    print(f"[接口] 蛋白 {prot.n_atoms} 原子 | 配体 {ligr} {lig.n_atoms} 原子 | "
          f"共 {len(u.trajectory)} 帧")

    p_elem, p_chg, p_name, p_rname, p_rid, p_chain = build_atom_arrays(u, prot)
    l_elem, l_chg, l_name, l_rname, l_rid, l_chain = build_atom_arrays(u, lig)

    hmap = attached_h_map(u)
    p_has_h = np.array([len(hmap.get(int(i), [])) > 0 for i in prot.indices])
    l_has_h = np.array([len(hmap.get(int(i), [])) > 0 for i in lig.indices])

    events = []
    rin = Counter()
    residue_contacts = Counter()
    picks = list(range(0, len(u.trajectory), max(1, args.stride)))
    if args.frames and len(picks) > args.frames:
        step = max(1, len(picks) // args.frames)
        picks = picks[::step][:args.frames]
    for fi, ts_idx in enumerate(picks):
        ts = u.trajectory[ts_idx]
        P = prot.positions
        L = lig.positions
        if cKDTree is not None:
            pairs = cKDTree(L).query_ball_point(P, args.cutoff)
        else:
            dd = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=-1)
            pairs = [np.where(dd[i] <= args.cutoff)[0] for i in range(len(P))]

        for i, js in enumerate(pairs):
            if len(js) == 0:
                continue
            for j in js:
                j = int(j)
                d = float(np.linalg.norm(P[i] - L[j]))
                if d > args.cutoff:
                    continue
                pe = p_elem[i]
                pc = float(p_chg[i])
                le = l_elem[j]
                lc = float(l_chg[j])
                prn = p_rname[i]
                lrn = l_rname[i] if False else l_rname[j]
                lrn = l_rname[j]

                cls = "Other Contact"
                status = "proximity-only"
                source = "no specific evidence"
                angle = None
                donor = ""
                acceptor = ""

                # --- 氢键: 需要 D-H...A 完整几何证据 ---
                hb = False
                if is_donor(pe, p_has_h[i]) and is_acceptor(le, lc, l_has_h[j]):
                    for hidx in hmap.get(int(prot.indices[i]), []):
                        hp = u.atoms[hidx].position
                        if float(np.linalg.norm(hp - L[j])) <= args.hbond:
                            ang = angle_deg(P[i], hp, L[j])
                            if ang is not None and ang >= args.hbond_angle:
                                hb = True
                                angle = ang
                                donor = p_name[i]
                                acceptor = l_name[j]
                                break
                if not hb and is_donor(le, l_has_h[j]) and is_acceptor(pe, pc, p_has_h[i]):
                    for hidx in hmap.get(int(lig.indices[j]), []):
                        hp = u.atoms[hidx].position
                        if float(np.linalg.norm(hp - P[i])) <= args.hbond:
                            ang = angle_deg(L[j], hp, P[i])
                            if ang is not None and ang >= args.hbond_angle:
                                hb = True
                                angle = ang
                                donor = l_name[j]
                                acceptor = p_name[i]
                                break
                if hb:
                    cls = "Hydrogen Bond"
                    status = "confirmed"
                    source = "D-H...A geometry"
                elif d <= args.salt and salt_bridge(prn, lrn, pc, lc):
                    cls = "Salt Bridge"
                    status = "confirmed"
                    source = "residue-charge + distance"
                elif d <= args.metal and pe in METAL_ELEMENTS:
                    cls = "Metal Coordination"
                    status = "confirmed"
                    source = "metal distance"
                elif pe in ("C", "S") and le in ("C", "S"):
                    cls = "Hydrophobic Contact"
                    status = "geometric"
                    source = "nonpolar heavy atoms"
                elif d <= args.salt and (
                    (prn in POS_RES and lc <= -0.4) or (lrn in POS_RES and pc <= -0.4)
                ):
                    cls = "Salt Bridge"
                    status = "geometric"
                    source = "residue identity + distance"

                events.append({
                    "Frame": ts_idx,
                    "Time_ns": round(float(ts.time) / 1000.0, 4),
                    "InteractionClass": cls,
                    "ValidationStatus": status,
                    "ClassificationSource": source,
                    "LigandResidue": f"{lrn}{int(l_rid[j])}",
                    "LigandAtom": l_name[j],
                    "LigandElement": le,
                    "LigandCharge": round(lc, 3),
                    "ProteinResidue": f"{prn}{int(p_rid[i])}",
                    "ProteinAtom": p_name[i],
                    "ProteinElement": pe,
                    "ProteinCharge": round(pc, 3),
                    "Distance": round(d, 3),
                    "Angle": round(angle, 1) if angle is not None else "",
                    "DonorAtom": donor,
                    "AcceptorAtom": acceptor,
                })
                key = (prn, int(p_rid[i]), p_chain[i], lrn, int(l_rid[j]))
                residue_contacts[key] += 1
                rin[(f"{p_chain[i]}:{prn}{int(p_rid[i])}",
                     f"{lrn}{int(l_rid[j])}")] += 1
        if (fi + 1) % 20 == 0 or fi == len(picks) - 1:
            print(f"  已处理 {fi+1}/{len(picks)} 帧, 事件 {len(events)} 条")

    n_frames = len(picks)
    ev_path = outdir / "interface_events.csv"
    with open(ev_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=EVENT_COLUMNS)
        w.writeheader()
        w.writerows(events)

    res_path = outdir / "interface_residues.csv"
    with open(res_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ProteinResidue", "ProteinChain", "LigandResidue",
                    "ContactCount", "ContactsPerFrame"])
        for (pn, pid, pc, ln, lid), c in residue_contacts.most_common():
            w.writerow([f"{pn}{pid}", pc, f"{ln}{lid}", c,
                        round(c / max(1, n_frames), 2)])

    rin_path = outdir / "interface_rin.csv"
    with open(rin_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ProteinResidue", "LigandResidue", "ContactCount",
                    "ContactsPerFrame", "FractionOfFrames"])
        for (a, b), c in rin.most_common():
            # 该残基对在多少比例的帧里出现过(每帧至多计 1 次, 这里用计数近似)
            w.writerow([a, b, c, round(c / max(1, n_frames), 2),
                        round(min(1.0, c / max(1, n_frames)), 4)])

    by_class = Counter(e["InteractionClass"] for e in events)
    by_status = Counter(e["ValidationStatus"] for e in events)
    summary = {
        "frames_analysed": n_frames,
        "total_events": len(events),
        "events_per_frame": round(len(events) / max(1, n_frames), 2),
        "ligand_resname": ligr,
        "protein_atoms": int(prot.n_atoms),
        "ligand_atoms": int(lig.n_atoms),
        "by_interaction_class": dict(by_class.most_common()),
        "by_validation_status": dict(by_status.most_common()),
        "interface_residues": len(residue_contacts),
        "top_interface_residues": [
            {"residue": f"{pn}{pid}:{pc}", "ligand": f"{ln}{lid}", "contacts": c}
            for (pn, pid, pc, ln, lid), c in residue_contacts.most_common(20)],
        "cutoffs": {"general": args.cutoff, "hbond": args.hbond,
                    "hbond_angle": args.hbond_angle, "salt_bridge": args.salt,
                    "metal": args.metal},
        "note": ("相互作用按拓扑信息(元素/电荷/残基名)与几何证据保守分类; "
                 "判不出具体类型的一律记为 Other Contact, 不臆测化学类型。"
                 "配体侧未做芳香环感知, 因此未单列 pi-pi 类别。"),
    }
    (outdir / "interface_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    plot_map(residue_contacts, outdir / "interface_map.png", n_frames)

    print()
    print(f"[接口] 事件 {len(events)} 条 / {n_frames} 帧 "
          f"(平均 {summary['events_per_frame']} 条每帧)")
    print("[接口] 分类统计:")
    for k, v in by_class.most_common():
        print(f"    {k:<24} {v}")
    print(f"[接口] 界面残基 {len(residue_contacts)} 个")
    print(f"[接口] 输出目录: {outdir}")
    return 0


def _find_cjk_font():
    """找一个能显示中文的 matplotlib 字体; 找不到返回 None(就用英文标签)。"""
    try:
        from matplotlib import font_manage
    except Exception:
        return None
    prefer = ["Noto Sans CJK SC", "Noto Sans CJK JP", "Noto Sans CJK TC",
              "WenQuanYi Zen Hei", "WenQuanYi Micro Hei", "Source Han Sans SC",
              "AR PL UMing CN", "Droid Sans Fallback"]
    have = {f.name for f in font_manager.fontManager.ttflist}
    for nm in prefer:
        if nm in have:
            return nm
    return None


def plot_map(residue_contacts, png, n_frames):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[接口] 画图依赖缺失, 跳过热图: {exc}", file=sys.stderr)
        return
    if not residue_contacts:
        return
    labels = []
    ligs = []
    for (pn, pid, pc, ln, lid) in residue_contacts:
        key = f"{pn}{pid}"
        if key not in labels:
            labels.append(key)
        lk = f"{ln}{lid}"
        if lk not in ligs:
            ligs.append(lk)
    labels = labels[:30]
    ligs = ligs[:20]
    lindex = {x: i for i, x in enumerate(ligs)}
    mat = np.zeros((len(labels), len(ligs)))
    for (pn, pid, pc, ln, lid), c in residue_contacts.items():
        key = f"{pn}{pid}"
        lk = f"{ln}{lid}"
        if key in labels and lk in lindex:
            mat[labels.index(key), lindex[lk]] = c
    if mat.sum() == 0:
        return
    fig, ax = plt.subplots(figsize=(max(6, len(ligs) * 0.4),
                                    max(4, len(labels) * 0.32)))
    im = ax.imshow(mat, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(ligs)))
    ax.set_xticklabels(ligs, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    # 图上用英文标签: matplotlib 默认字体(DejaVu Sans)不含 CJK,
    # 中文会变成方块; 有 CJK 字体时再切回去
    cn_font = _find_cjk_font()
    if cn_font:
        matplotlib.rcParams["font.sans-serif"] = [cn_font]
        matplotlib.rcParams["axes.unicode_minus"] = False
        xlabel, ylabel, title = "Ligand residue", "Protein interface residues", \
            f"Interface contact frequency ({n_frames} frames)"
    else:
        xlabel, ylabel, title = "Ligand residue", "Protein interface residues", \
            f"Interface contact frequency ({n_frames} frames)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="contacts")
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print(f"[接口] 热图 {png}")


if __name__ == "__main__":
    sys.exit(main())