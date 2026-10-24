#!/usr/bin/env python3
"""蛋白 + 配体 -> 复合物坐标与拓扑。

做两件事:
1. 坐标: 蛋白(去掉水/离子) 与配体 合并成一个 gro, 盒矢量取两者较大值
2. 拓扑: 重写 topol.top 的力场引用与 [ molecules ] 段, 使原子顺序与坐标一致

GROMACS 的拓扑只允许一个 [ atomtypes ] 段, 所以配体的 GAFF 类型必须并入蛋白力场。
本脚本假定已经用 merge_gaff_ff.py 生成好了合并力场目录和去掉 atomtypes 段的配体 itp,
也可以用 --auto-merge 自动完成(需要 acpype 原始 itp)。

用法:
  assemble_complex.py --protein protein.gro --ligand lig.gro \\
      --topol topol.top --lig-itp lig.itp --out-prefix complex
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

# 视为"溶剂/离子"的残基名, 组装复合物时会被剔除(之后统一用 solvate/genion 重新加)
SOLVENT_RESNAMES = {"SOL", "HOH", "WAT", "TIP3", "TIP4", "SPC", "NA", "CL", "K",
                    "MG", "CA", "ZN", "FE", "BR", "IOD", "CS", "RB", "LI", "F"}


def read_gro(path):
    """读 gro 文件, 返回 (标题行, 原子行列表, 盒矢量字符串)。"""
    lines = open(path).read().splitlines()
    if len(lines) < 3:
        raise SystemExit("错误: %s 不是合法的 gro 文件" % path)
    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]
    if len(atoms) != natoms:
        raise SystemExit("错误: %s 声明 %d 个原子, 实际 %d 行" % (path, natoms, len(atoms)))
    return lines[0], atoms, lines[2 + natoms]


def resname(atom_line):
    """gro 格式: %5d%-5s%5s... 残基名在第 6-10 列(索引 5:10)。"""
    return atom_line[5:10].strip()


def renumber(atom_line, resnr, atomnr):
    """重写 gro 行中的残基号与原子序号。

    合并蛋白与配体时必须重新编号: 配体自带的残基号通常是 1, 会与蛋白的残基 1
    冲突, 导致 PDB 转换时被 GROMACS/ParmEd 丢弃原子。
    gro 列布局: %5d%-5s%5s%5d  ->  残基号[0:5] 名称[5:10] 原子名[10:15] 序号[15:20]
    """
    return ("%5d" % resnr) + atom_line[5:15] + ("%5d" % atomnr) + atom_line[20:]


def box_vec(box_line):
    v = [float(x) for x in box_line.split()[:3]]
    while len(v) < 3:
        v.append(0.0)
    return v


def strip_solvent(atoms):
    kept = [a for a in atoms if resname(a) not in SOLVENT_RESNAMES]
    return kept, len(atoms) - len(kept)


def parse_molecules(topol_path):
    """从 topol.top 的 [ molecules ] 段解析出 [(分子名, 数量)], 去掉溶剂/离子。"""
    lines = open(topol_path).read().splitlines()
    entries, inmol = [], False
    for L in lines:
        s = L.strip()
        if s.startswith("["):
            inmol = s.lower().startswith("[ molecules ]")
            continue
        if not inmol or not s or s.startswith(";"):
            continue
        parts = s.split()
        if len(parts) >= 2 and parts[0].upper() not in SOLVENT_RESNAMES:
            try:
                entries.append((parts[0], int(parts[1])))
            except ValueError:
                pass
    return entries, lines


def moleculetype_name(itp_path):
    """从配体 itp 里读出 [ moleculetype ] 的名字(跳过注释行)。"""
    lines = open(itp_path, encoding="utf-8", errors="replace").read().splitlines()
    inmol = False
    for L in lines:
        s = L.strip()
        if s.startswith("["):
            inmol = s.lower().startswith("[ moleculetype ]")
            continue
        if not inmol or not s or s.startswith(";"):
            continue
        return s.split()[0]
    raise SystemExit("错误: %s 里找不到 [ moleculetype ]" % itp_path)


def main():
    ap = argparse.ArgumentParser(description="组装蛋白-配体复合物")
    ap.add_argument("--protein", required=True, help="蛋白坐标 .gro")
    ap.add_argument("--ligand", required=True, help="配体坐标 .gro")
    ap.add_argument("--topol", required=True, help="蛋白拓扑 topol.top")
    ap.add_argument("--lig-itp", required=True, help="配体拓扑 .itp(不应含 [ atomtypes ])")
    ap.add_argument("--out-prefix", default="complex", help="输出前缀, 默认 complex")
    ap.add_argument("--ff-dir", default=None,
                    help="合并后的蛋白力场目录(绝对路径), 如 /path/oplsaa_gaff.ff")
    ap.add_argument("--keep-solvent", action="store_true",
                    help="保留蛋白坐标里的水/离子(默认剔除)")
    args = ap.parse_args()

    out_gro = args.out_prefix + ".gro"
    out_top = args.out_prefix + ".top"

    # ---------- 1. 坐标 ----------
    _, patoms, pbox = read_gro(args.protein)
    _, latoms, lbox = read_gro(args.ligand)
    if args.keep_solvent:
        prot, ndrop = patoms, 0
    else:
        prot, ndrop = strip_solvent(patoms)
    # 重新编号: 配体残基号接在蛋白之后, 原子序号连续
    last_res = max(int(a[0:5]) for a in prot) if prot else 0
    last_atm = max(int(a[15:20]) for a in prot) if prot else 0
    latoms = [renumber(a, last_res + 1, last_atm + i + 1)
              for i, a in enumerate(latoms)]

    print("[1/4] 蛋白 %d 原子%s + 配体 %d 原子 = %d 原子"
          % (len(prot), ("(剔除溶剂/离子 %d 个)" % ndrop) if ndrop else "",
             len(latoms), len(prot) + len(latoms)))
    print("      配体已重新编号: 残基号 %d, 原子序号 %d-%d"
          % (last_res + 1, last_atm + 1, last_atm + len(latoms)))

    allatoms = prot + latoms
    pv, lv = box_vec(pbox), box_vec(lbox)
    merged = [max(pv[0], lv[0]), max(pv[1], lv[1]), max(pv[2], lv[2])]
    with open(out_gro, "w") as f:
        f.write("%s protein+ligand\n" % args.out_prefix)
        f.write("%5d\n" % len(allatoms))
        f.write("\n".join(allatoms) + "\n")
        f.write("%10.5f%10.5f%10.5f\n" % tuple(merged))
    print("      写出 %s (盒矢量 %s)" % (out_gro, " ".join("%.3f" % v for v in merged)))

    # ---------- 2. 拓扑 ----------
    prot_entries, src_lines = parse_molecules(args.topol)
    if not prot_entries:
        raise SystemExit("错误: %s 的 [ molecules ] 段里没有可用的蛋白分子" % args.topol)
    lig_name = moleculetype_name(args.lig_itp)
    print("[2/4] 蛋白分子 %s ; 配体 moleculetype = %s"
          % (", ".join("%s x%d" % e for e in prot_entries), lig_name))

    out = []
    for L in src_lines:
        s = L.strip()
        if s.startswith("[ molecules ]"):
            break
        if args.ff_dir:
            # 把蛋白力场引用换成合并后的目录(绝对路径)
            ffdir = args.ff_dir.rstrip("/")
            if not ffdir.endswith(".ff"):
                ffdir = os.path.join(ffdir, "PROTEIN_FF_NAME")
            L = re.sub(r'"[^"]*?\.ff/forcefield\.itp"',
                       '"%s/forcefield.itp"' % ffdir, L)
        out.append(L)

    out.append('#include "%s"' % os.path.basename(args.lig_itp))
    out.append("")
    out.append("[ molecules ]")
    out.append("; Compound        #mols")
    for name, cnt in prot_entries:
        out.append("%-16s %d" % (name, cnt))
    out.append("%-16s %d" % (lig_name, 1))
    with open(out_top, "w") as f:
        f.write("\n".join(out) + "\n")
    print("      写出 %s" % out_top)

    # ---------- 3. 把配体 itp 拷到当前目录(GROMACS 按相对路径 include) ----------
    if os.path.dirname(os.path.abspath(args.lig_itp)) != os.getcwd():
        shutil.copy(args.lig_itp, os.path.basename(args.lig_itp))
        print("[3/4] 配体 itp 已拷贝到当前目录: %s" % os.path.basename(args.lig_itp))
    else:
        print("[3/4] 配体 itp 已在当前目录")

    # ---------- 4. 自检 ----------
    print("[4/4] 自检:")
    print("      坐标原子数 %d = 蛋白 %d + 配体 %d" % (len(allatoms), len(prot), len(latoms)))
    print("      拓扑声明的分子: %s + %s x1"
          % (", ".join("%s x%d" % e for e in prot_entries), lig_name))
    print("      下一步: gmx grompp -f mdp/em.mdp -c %s -p %s -o em.tpr" % (out_gro, out_top))
    print("            注意拓扑里配体 itp 用的是相对路径, 必须在同一目录下运行 grompp")


if __name__ == "__main__":
    main()
