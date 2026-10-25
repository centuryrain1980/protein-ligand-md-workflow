#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""protein-ligand-md-workflow 统一命令行入口。

把 skill 里的各个 shell 脚本串成一个可被 Agent / 人直接调用的流水线。

用法:
    mdflow.py <命令> [选项]
    mdflow.py info          # 打印流水线说明
    mdflow.py --help        # 完整帮助

命令:
    env        环境自检(逐项打印 OK/WARN/FAIL)
    install    一键安装全部运行环境
    verify     校验当前体系是否具备跑 MD 的条件
    dock       分子对接(gnina): 配体进口袋
    prepare    蛋白 + 配体 -> 复合物(含 GAFF 力场合并)
    solvate    溶剂化 + 加离子
    simulate   EM -> NVT -> NPT -> 生产 MD
    analyze    轨迹分析(出图)
    mmpbsa     MM/PBSA 结合自由能 + 能量分解
    report     汇总成 PDF 报告
    all        按顺序跑完整流水线(verify->prepare->solvate->simulate->analyze->report)

通用选项:
    --workdir DIR     工作目录(默认当前目录)
    --dry-run         只打印将要执行的命令, 不真的执行
    -v/--verbose      打印子命令的完整输出

示例:
    mdflow.py prepare --receptor rec.pdb --ligand-file lig.mol2 --workdir proj
    mdflow.py solvate --workdir proj --box 3.0
    mdflow.py simulate --workdir proj --nsteps 50000000 --resume
    mdflow.py report  --workdir proj --title "KRAS G12D - Sotorasib"
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPTS = SKILL_DIR / "scripts"

GREEN, YELLOW, RED, BLUE, BOLD, NC = "\033[32m", "\033[33m", "\033[31m", "\033[34m", "\033[1m", "\033[0m"
if not sys.stdout.isatty():
    GREEN = YELLOW = RED = BLUE = BOLD = NC = ""


def log(msg):
    print(f"{BLUE}[mdflow]{NC} {msg}")


def ok(msg):
    print(f"  {GREEN}[OK]{NC}   {msg}")


def warn(msg):
    print(f"  {YELLOW}[WARN]{NC} {msg}")


def bad(msg):
    print(f"  {RED}[FAIL]{NC} {msg}")


def die(msg, code=1):
    print(f"  {RED}[FAIL]{NC} {msg}", file=sys.stderr)
    sys.exit(code)


def which(name):
    return shutil.which(name)


def run(cmd, dry=False, verbose=False, shell=False, cwd=None, env=None):
    """执行命令。返回 (rc, 输出文本)。"""
    printable = cmd if isinstance(cmd, str) else " ".join(cmd)
    if dry:
        print(f"  {YELLOW}[dry-run]{NC} {printable}")
        return 0, ""
    log(f"→ {printable}")
    e = dict(os.environ)
    if env:
        e.update(env)
    try:
        p = subprocess.run(cmd, shell=shell, cwd=cwd, env=e, text=True,
                           capture_output=True, errors="replace")
    except FileNotFoundError:
        die(f"找不到命令: {cmd}")
    out = (p.stdout or "") + (p.stderr or "")
    if verbose and out.strip():
        print(out.rstrip())
    if p.returncode != 0:
        if not verbose:
            tail = "\n".join(out.strip().splitlines()[-12:])
            if tail:
                print(f"  {RED}--- 输出尾部 ---{NC}\n{tail}")
        die(f"命令失败 (exit={p.returncode}): {printable}", p.returncode)
    return p.returncode, out


# ----------------------------------------------------------------------
def cmd_env(a):
    run([str(SCRIPTS / "setup_env.sh")], a.dry_run, a.verbose)
    return 0


def cmd_install(a):
    run(["bash", str(SCRIPTS / "install_all.sh")], a.dry_run, a.verbose)
    return 0


def cmd_verify(a):
    w = Path(a.workdir).resolve()
    log(f"校验体系: {w}")
    need = [("蛋白拓扑", "topol.top"), ("蛋白坐标", "protein.gro"),
            ("复合物坐标", "complex.gro"), ("复合物拓扑", "complex.top"),
            ("溶剂化结构", "solv.gro"), ("加离子结构", "solv_ions.gro"),
            ("生产拓扑", "md.tpr"), ("生产轨迹", "md.xtc")]
    missing = []
    for label, fn in need:
        p = w / fn
        if p.is_file():
            ok(f"{label}: {fn} ({p.stat().st_size} bytes)")
        else:
            miss = {"topol.top": "跑 prepare", "protein.gro": "跑 pdb2gmx",
                    "complex.gro": "跑 prepare", "complex.top": "跑 prepare",
                    "solv.gro": "跑 solvate", "solv_ions.gro": "跑 solvate",
                    "md.tpr": "跑 simulate", "md.xtc": "跑 simulate"}
            warn(f"{label}: 缺 {fn} → {miss.get(fn, '?')}")
            missing.append(fn)
    # 工具
    for t in ("gmx", "dit", "obabel"):
        (ok if which(t) else warn)(f"工具 {t}: {which(t) or '未安装'}")
    if a.ff:
        ffd = Path(os.environ.get("GMXDATA", "/usr/local/share/gromacs/top")) / f"{a.ff}.ff"
        (ok if ffd.is_dir() else warn)(f"力场 {a.ff}: {ffd if ffd.is_dir() else '不存在'}")
    if missing:
        print(f"\n  结论: 体系还不完整, 缺 {len(missing)} 项")
        return 1
    print("\n  结论: 体系完整, 可以跑 simulate")
    return 0


def cmd_dock(a):
    if not a.receptor or not a.ligand_file:
        die("dock 需要 --receptor 和 --ligand-file")
    cmd = ["bash", str(SCRIPTS / "dock_ligand.sh"), a.receptor, a.ligand_file]
    if a.center:
        cmd += [str(x) for x in a.center]
        if a.box:
            cmd.append(str(a.box))
    run(cmd, a.dry_run, a.verbose, cwd=a.workdir)
    return 0


def _ff_candidates(a):
    """解析力场候选列表。

    --ff-list 语法: "amber99sb:tip3p,charmm36:tip3p,charmm27:tips3p"
    每项是 力场[:水模型]; 不写水模型时沿用 --water。
    为空则退回单一候选 (--ff, --water)。
    """
    if not a.ff_list:
        return [(a.ff, a.water)]
    out = []
    for item in a.ff_list.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            ff, wat = item.split(":", 1)
            out.append((ff.strip(), wat.strip() or a.water))
        else:
            out.append((item, a.water))
    return out or [(a.ff, a.water)]


def cmd_prepare(a):
    if not a.receptor or not a.ligand_file:
        die("prepare 需要 --receptor 和 --ligand-file")
    w = Path(a.workdir).resolve()
    # 先 pdb2gmx(力场候选逐个尝试, 前一个失败就换下一个)
    topol = w / "topol.top"
    cands = _ff_candidates(a)
    if not topol.is_file() or a.force:
        last = ""
        for idx, (ff, wat) in enumerate(cands):
            log(f"pdb2gmx 建蛋白拓扑 (候选 {idx+1}/{len(cands)}: {ff} + {wat})")
            cmd = ["gmx", "pdb2gmx", "-f", a.receptor, "-o", "protein.gro",
                   "-i", "posre.itp", "-p", "topol.top",
                   "-ff", ff, "-water", wat, "-ignh"]
            if a.dry_run:
                run(cmd, True, a.verbose, cwd=w)
                break
            p = subprocess.run(cmd, cwd=w, text=True, capture_output=True,
                               errors="replace")
            out = (p.stdout or "") + (p.stderr or "")
            if p.returncode == 0:
                ok(f"力场 {ff} + 水模型 {wat} 可用")
                a.ff, a.water = ff, wat
                break
            last = out
            warn(f"力场 {ff} + 水模型 {wat} 失败, 尝试下一个候选")
            if a.verbose:
                print("\n".join(out.strip().splitlines()[-8:]))
        else:
            bad("所有力场候选都失败, 最后一次的报错:")
            print("\n".join(last.strip().splitlines()[-12:]))
            return 1
    else:
        log("topol.top 已存在, 跳过 pdb2gmx (加 --force 可重建)")
    # 再 setup_ligand.sh
    env = {"PROT_FF": a.ff, "PROT_WATER": a.water}
    cmd = ["bash", str(SCRIPTS / "setup_ligand.sh"), a.ligand_file,
           "topol.top", "protein.gro", "complex"]
    run(cmd, a.dry_run, a.verbose, cwd=w, env=env)
    return 0


def cmd_solvate(a):
    w = Path(a.workdir).resolve()
    if not (w / "complex.gro").is_file():
        die("找不到 complex.gro, 先跑 prepare")
    box = str(a.box)
    run(["gmx", "editconf", "-f", "complex.gro", "-o", "complex_c.gro",
         "-c", "-box", box, box, box], a.dry_run, a.verbose, cwd=w)
    run(["gmx", "solvate", "-cp", "complex_c.gro",
         "-cs", "/usr/local/share/gromacs/top/spc216.gro",
         "-box", box, box, box, "-o", "solv.gro", "-p", "complex.top"],
        a.dry_run, a.verbose, cwd=w)
    run(["gmx", "grompp", "-f", str(SKILL_DIR / "assets/mdp/ions.mdp"),
         "-c", "solv.gro", "-p", "complex.top", "-o", "ions.tpr"],
        a.dry_run, a.verbose, cwd=w)
    # genion 需要交互喂组名
    if a.dry_run:
        log("[dry-run] printf 'SOL\\n' | gmx genion -s ions.tpr -o solv_ions.gro -p complex.top "
            "-pname NA -nname CL -neutral")
    else:
        log("→ printf 'SOL\\n' | gmx genion ...")
        p = subprocess.run("printf 'SOL\\n' | gmx genion -s ions.tpr -o solv_ions.gro "
                           "-p complex.top -pname NA -nname CL -neutral",
                           shell=True, cwd=w, text=True, capture_output=True,
                           errors="replace")
        if p.returncode != 0:
            print((p.stdout or "") + (p.stderr or ""))
            die("genion 失败")
        ok("solv_ions.gro 已生成")
    return 0


def cmd_simulate(a):
    w = Path(a.workdir).resolve()
    cmd = ["bash", str(SCRIPTS / "run_simulation.sh"), str(w), str(a.nsteps)]
    env = {"NTOMP": str(a.ntomp)}
    if a.resume:
        cmd.append("--resume")
    if a.plan:
        cmd += ["--plan", a.plan]
    run(cmd, a.dry_run, a.verbose, cwd=w, env=env)
    return 0


def cmd_analyze(a):
    w = Path(a.workdir).resolve()
    cmd = ["bash", str(SCRIPTS / "analyze.sh"), str(w), a.prot, a.lig]
    run(cmd, a.dry_run, a.verbose, cwd=w)
    return 0


def cmd_mmpbsa(a):
    w = Path(a.workdir).resolve()
    cmd = ["bash", str(SCRIPTS / "run_mmpbsa.sh"), str(w), a.prot, a.lig,
           str(a.startframe), str(a.endframe)]
    run(cmd, a.dry_run, a.verbose, cwd=w)
    return 0


def cmd_report(a):
    w = Path(a.workdir).resolve()
    cmd = ["python3", str(SCRIPTS / "make_report.py"), str(w)]
    if a.title:
        cmd += ["--title", a.title]
    if a.ligand_file:
        cmd += ["--ligand-file", a.ligand_file]
    if a.html:
        cmd.append("--html")
    run(cmd, a.dry_run, a.verbose, cwd=w)
    return 0


def cmd_all(a):
    steps = [
        ("校验体系", cmd_verify),
        ("建体系(蛋白+配体)", cmd_prepare),
        ("溶剂化+加离子", cmd_solvate),
        ("平衡+生产 MD", cmd_simulate),
        ("轨迹分析", cmd_analyze),
        ("生成报告", cmd_report),
    ]
    for name, fn in steps:
        print(f"\n{BOLD}===== {name} ====={NC}")
        rc = fn(a)
        if rc not in (0, None):
            print(f"  {YELLOW}[SKIP]{NC} 步骤 '{name}' 未完成, 后续步骤可能失败")
    return 0


def cmd_info(a):
    print(f"""{BOLD}protein-ligand-md-workflow 流水线{NC}

  verify    校验体系完整性
  prepare   蛋白(pdb2gmx) + 配体(GAFF) -> 复合物
  solvate   建盒、加水、加离子
  simulate  EM -> NVT -> NPT -> 生产 MD
  analyze   RMSD/RMSF/Rg/SASA/Rama/二级结构/能量, 出图
  mmpbsa    MM/PBSA 结合自由能与逐残基分解
  report    汇总成 PDF 报告 (analysis/report.pdf)

典型用法:
  mdflow.py prepare --receptor rec.pdb --ligand-file lig.mol2 --workdir proj
  mdflow.py solvate --workdir proj --box 3.0
  mdflow.py simulate --workdir proj --nsteps 50000000 --resume
  mdflow.py analyze --workdir proj --prot Protein --lig LIG
  mdflow.py report  --workdir proj --title "我的体系"

提示: 任何命令加 --dry-run 可以先看会执行什么, 不真正跑。
""")
    return 0


COMMANDS = {
    "env": cmd_env, "install": cmd_install, "verify": cmd_verify,
    "dock": cmd_dock, "prepare": cmd_prepare, "solvate": cmd_solvate,
    "simulate": cmd_simulate, "analyze": cmd_analyze, "mmpbsa": cmd_mmpbsa,
    "report": cmd_report, "all": cmd_all, "info": cmd_info,
}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="protein-ligand-md-workflow 统一命令行入口",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("command", choices=sorted(COMMANDS), help="要执行的步骤")
    ap.add_argument("--workdir", default=".", help="工作目录(默认当前目录)")
    ap.add_argument("--receptor", help="受体 PDB")
    ap.add_argument("--ligand-file", help="配体 mol2/sdf")
    ap.add_argument("--ligand", help="配体文件(同 --ligand-file)")
    ap.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                    help="对接盒子中心")
    ap.add_argument("--box", type=float, default=3.0, help="溶剂化盒子边长 nm")
    ap.add_argument("--ff", default="amber99sb", help="蛋白力场(默认 amber99sb)")
    ap.add_argument("--water", default="tip3p", help="水模型(默认 tip3p)")
    ap.add_argument("--ff-list", default=None,
                    help="力场候选列表(逗号分隔), 如 'amber99sb:tip3p,charmm36:tip3p'; "
                         "前一个建拓扑失败就自动换下一个")
    ap.add_argument("--prot", default="Protein", help="蛋白组名")
    ap.add_argument("--lig", default="LIG", help="配体组名")
    ap.add_argument("--nsteps", type=int, default=50000000, help="生产 MD 步数")
    ap.add_argument("--ntomp", type=int, default=6, help="OpenMP 线程数")
    ap.add_argument("--resume", action="store_true", help="从检查点续跑生产 MD")
    ap.add_argument("--plan", default=None,
                    help="JSON 多阶段计划文件(见 assets/plan/*.json), 会覆盖默认四阶段")
    ap.add_argument("--force", action="store_true", help="强制重建(即使已存在)")
    ap.add_argument("--startframe", type=int, default=1, help="MM/PBSA 起始帧")
    ap.add_argument("--endframe", type=int, default=50, help="MM/PBSA 结束帧")
    ap.add_argument("--title", help="报告标题")
    ap.add_argument("--html", action="store_true", help="报告输出 HTML 而非 PDF")
    ap.add_argument("--dry-run", action="store_true", help="只打印命令不执行")
    ap.add_argument("-v", "--verbose", action="store_true", help="打印子命令完整输出")
    a = ap.parse_args(argv)

    if not a.ligand_file and a.ligand:
        a.ligand_file = a.ligand
    a.workdir = Path(a.workdir).expanduser().resolve()
    if a.command not in ("info", "env", "install") and not a.workdir.is_dir():
        die(f"工作目录不存在: {a.workdir}")
    return COMMANDS[a.command](a) or 0


if __name__ == "__main__":
    sys.exit(main())