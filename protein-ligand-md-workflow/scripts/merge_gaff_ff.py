#!/usr/bin/env python3
"""把 GAFF 配体原子类型并入蛋白力场目录, 供 GROMACS 2025 使用。

背景: GROMACS 的拓扑只允许一个 [ atomtypes ] 段, 而 acpype 生成的配体 itp 自带
一段, 直接 include 会报 "Syntax error"。GROMACS 2025 的力场目录把原子类型分散在
两个文件里, 必须同时写入, 否则会报 "Atomtype XXX not found":
  - ffnonbonded.itp : [ atomtypes ] 段, 含 bond_type/charge/ptype/sigma/epsilon
  - atomtypes.atp   : name + mass 列表

用法:
  merge_gaff_ff.py --lig-itp lig_GMX.itp --ff-src /usr/local/share/gromacs/top/oplsaa.ff \
                    --out-dir ./oplsaa_gaff.ff
  merge_gaff_ff.py --lig-itp lig_GMX.itp --ff-src ... --out-dir ... --strip-lig-itp lig.itp
"""
import argparse
import os
import re
import shutil
import sys


def read_gaff_types(lig_itp):
    """从 acpype 生成的 itp 里取出 GAFF 原子类型行。"""
    txt = open(lig_itp, encoding="utf-8", errors="replace").read()
    if "[ atomtypes ]" not in txt:
        raise SystemExit("错误: %s 里没有 [ atomtypes ] 段" % lig_itp)
    blk = txt.split("[ atomtypes ]", 1)[1].split("]", 1)[0]
    rows = []
    for line in blk.splitlines():
        s = line.strip()
        if not s or s.startswith(";") or s.startswith("["):
            continue
        p = s.split()
        if len(p) < 3:
            continue
        # acpype 格式: name bond_type mass charge ptype sigma epsilon [; amb]
        rows.append((p[0], p[2], s))
    if not rows:
        raise SystemExit("错误: 未从 %s 解析出任何原子类型" % lig_itp)
    return rows


def main():
    ap = argparse.ArgumentParser(description="把 GAFF 原子类型并入 GROMACS 力场目录")
    ap.add_argument("--lig-itp", required=True, help="acpype 生成的配体 itp")
    ap.add_argument("--ff-src", required=True, help="蛋白力场目录, 如 .../top/oplsaa.ff")
    ap.add_argument("--out-dir", required=True, help="输出的合并后力场目录")
    ap.add_argument("--strip-lig-itp", default=None,
                    help="同时生成去掉 [ atomtypes ] 段的配体 itp 到此路径")
    args = ap.parse_args()

    rows = read_gaff_types(args.lig_itp)
    print("[1/4] 解析到 %d 个 GAFF 原子类型: %s"
          % (len(rows), ", ".join(r[0] for r in rows)))

    # 2. 复制力场目录
    if os.path.exists(args.out_dir):
        shutil.rmtree(args.out_dir)
    shutil.copytree(args.ff_src, args.out_dir)
    print("[2/4] 已复制 %s -> %s" % (args.ff_src, args.out_dir))

    # 3. 写入原子类型
    fn_nb = os.path.join(args.out_dir, "ffnonbonded.itp")
    with open(fn_nb, "a", encoding="utf-8") as f:
        f.write("\n; ---- GAFF atom types for ligand (from acpype) ----\n")
        for _, _, s in rows:
            f.write(s + "\n")
    print("[3/4] 已追加到 %s" % fn_nb)

    fn_atp = os.path.join(args.out_dir, "atomtypes.atp")
    if os.path.exists(fn_atp):
        with open(fn_atp, "a", encoding="utf-8") as f:
            for name, mass, _ in rows:
                f.write(" %-10s %s  ;\n" % (name, mass))
        print("      已追加到 %s" % fn_atp)
    else:
        print("      警告: %s 不存在, 跳过(某些力场没有该文件)" % fn_atp)

    # 4. 生成去掉 [ atomtypes ] 的配体 itp
    if args.strip_lig_itp:
        txt = open(args.lig_itp, encoding="utf-8", errors="replace").read()
        out = re.sub(r"\n[ \t]*\[ atomtypes \].*?(?=\n[ \t]*\[)", "\n", txt, flags=re.S)
        open(args.strip_lig_itp, "w", encoding="utf-8").write(out)
        print("[4/4] 已生成 %s (原子类型段已移除: %d -> %d 个)"
              % (args.strip_lig_itp, txt.count("[ atomtypes ]"), out.count("[ atomtypes ]")))
    else:
        print("[4/4] 未指定 --strip-lig-itp, 跳过")


if __name__ == "__main__":
    main()
