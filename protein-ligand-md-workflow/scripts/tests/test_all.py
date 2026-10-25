#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""protein-ligand-md-workflow 回归测试。

用法:
    python3 scripts/tests/test_all.py            # 跑全部
    python3 scripts/tests/test_all.py --fast     # 跳过需要跑 GROMACS 的慢测试
    python3 scripts/tests/test_all.py -k mdp     # 只跑名字里含 "mdp" 的用例

退出码 0 = 全通过, 1 = 有失败。
每个用例都是独立断言, 失败会打印具体原因。
"""

from __future__ import annotations

import argparse
import io
import json
import os
import re
import subprocess
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
SKILL = TESTS_DIR.parent.parent

GREEN, RED, YELLOW, BOLD, NC = "\033[32m", "\033[31m", "\033[33m", "\033[1m", "\033[0m"
if not sys.stdout.isatty():
    GREEN = RED = YELLOW = BOLD = NC = ""

_RESULTS = []


def case(name):
    """装饰器: 注册一个测试用例。"""
    def deco(fn):
        fn._case_name = name
        return fn
    return deco


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


def read(path):
    return io.open(path, encoding="utf-8", errors="replace").read()


# =========================================================================
# 1. 脚本语法
# =========================================================================
@case("语法: 所有 shell 脚本通过 bash -n")
def t_bash_syntax():
    bad = []
    for f in sorted((SKILL / "scripts").glob("*.sh")):
        p = subprocess.run(["bash", "-n", str(f)], capture_output=True, text=True)
        if p.returncode != 0:
            bad.append(f"{f.name}: {p.stderr.strip()[:120]}")
    check(not bad, "shell 脚本语法错误: " + "; ".join(bad))


@case("导入: 所有 Python 脚本能真正 import(抓被换行/吞字符造成的静默损坏)")
def t_py_import():
    """py_compile 只能抓语法错; 像 'from collections import Counte' 这种
    被吞掉一个字符的情况语法是合法的, 只有真正 import 才会暴露。"""
    import importlib.util
    bad = []
    for f in sorted(list((SKILL / "scripts").glob("*.py"))
                    + list((SKILL / "scripts/tests").glob("*.py"))):
        try:
            spec = importlib.util.spec_from_file_location(f"_t_{f.stem}", str(f))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
        except Exception as exc:
            bad.append(f"{f.name}: {type(exc).__name__}: {exc}")
    check(not bad, "脚本 import 失败: " + "; ".join(bad))


@case("语法: 所有 Python 脚本可编译")
def t_py_syntax():
    bad = []
    for f in sorted((SKILL / "scripts").glob("*.py")):
        p = subprocess.run([sys.executable, "-m", "py_compile", str(f)],
                           capture_output=True, text=True)
        if p.returncode != 0:
            bad.append(f"{f.name}: {p.stderr.strip()[:160]}")
        try:
            os.remove(f.parent / "__pycache__" / (f.stem + ".cpython-*.pyc"))
        except Exception:
            pass
    check(not bad, "Python 脚本编译失败: " + "; ".join(bad))


# =========================================================================
# 2. mdp 模板的关键参数
# =========================================================================
@case("mdp: CHARMM 模板含官方强制的 10 项参数")
def t_charmm_mdp():
    # EM / ions 阶段不做积分, 不需要 constraints; 其余 9 项所有阶段都要
    md_only = {"constraints": "h-bonds"}
    need = {
        "cutoff-scheme": "Verlet",
        "vdwtype": "cutoff",
        "vdw-modifier": "force-switch",
        "rlist": "1.2",
        "rvdw": "1.2",
        "rvdw-switch": "1.0",
        "coulombtype": "PME",
        "rcoulomb": "1.2",
        "DispCorr": "no",
    }
    for f in sorted((SKILL / "assets/mdp").glob("charmm36_*.mdp")):
        txt = read(f)
        want = dict(need)
        if not f.name.startswith(("charmm36_em", "charmm36_ions")):
            want.update(md_only)
        for key, val in want.items():
            m = re.search(r"^\s*%s\s*=\s*(\S+)" % re.escape(key), txt, re.M)
            check(m is not None, f"{f.name} 缺少 {key}")
            check(m.group(1) == val,
                  f"{f.name} 的 {key} 应为 {val}, 实际 {m.group(1)}")


@case("mdp: 通用模板用 Linear+System 消除 GROMACS 2025 警告")
def t_generic_mdp():
    for f in sorted((SKILL / "assets/mdp").glob("*.mdp")):
        if f.name.startswith("charmm36_"):
            continue
        txt = read(f)
        if "comm-mode" in txt:
            m = re.search(r"^\s*comm-mode\s*=\s*(\S+)", txt, re.M)
            check(m and m.group(1) == "Linear",
                  f"{f.name}: comm-mode 应为 Linear(Angular 会触发 2025 警告)")
        if "comm-grps" in txt:
            m = re.search(r"^\s*comm-grps\s*=\s*(\S+)", txt, re.M)
            check(m and m.group(1) == "System",
                  f"{f.name}: comm-grps 应为 System")


@case("mdp: CHARMM 模板的 rcoulomb 与 rvdw 相等(Verlet+PME 要求)")
def t_charmm_cutoff():
    for f in sorted((SKILL / "assets/mdp").glob("charmm36_*.mdp")):
        txt = read(f)
        rc = re.search(r"^\s*rcoulomb\s*=\s*(\S+)", txt, re.M)
        rv = re.search(r"^\s*rvdw\s*=\s*(\S+)", txt, re.M)
        check(rc and rv, f"{f.name} 缺少 rcoulomb/rvdw")
        check(rc.group(1) == rv.group(1),
              f"{f.name}: rcoulomb({rc.group(1)}) != rvdw({rv.group(1)}) "
              f"会触发 'rcoulomb!=rvdw is not supported'")


# =========================================================================
# 3. MM/PBSA 输入模板
# =========================================================================
@case("mmpbsa: 生成的 mmpbsa.in 不含 legacy forcefields 变量")
def t_mmpbsa_no_legacy():
    txt = read(SKILL / "scripts/run_mmpbsa.sh")
    body = txt.split("cat > mmpbsa.in <<EOF", 1)
    check(len(body) == 2, "run_mmpbsa.sh 里找不到 mmpbsa.in 模板")
    tmpl = body[1].split("EOF", 1)[0]
    for var in ("forcefields", "protein_forcefield", "ligand_forcefield",
                "ions_parameters"):
        check(not re.search(r"^\s*%s\s*=" % var, tmpl, re.M),
              f"mmpbsa.in 模板里不应再有 legacy 变量 {var}")
    check(re.search(r"^\s*source_force_field\s*=", tmpl, re.M),
          "mmpbsa.in 模板里应有 source_force_field")


@case("mmpbsa: print_res 用官方推荐的 within 6")
def t_mmpbsa_print_res():
    txt = read(SKILL / "scripts/run_mmpbsa.sh")
    check("within 6" in txt, "run_mmpbsa.sh 里应有 PRINT_RES 默认值 within 6")
    check(re.search(r'PRINT_RES="\$\{PRINT_RES:-within 6\}"', txt),
          "PRINT_RES 默认值应为 within 6")


@case("mmpbsa: igb 与 PBRadii 配对正确(igb=8 -> PBRadii=4/mbondi3)")
def t_mmpbsa_igb():
    txt = read(SKILL / "scripts/run_mmpbsa.sh")
    check(re.search(r'IGB="\$\{IGB:-8\}"', txt), "IGB 默认应为 8 (GB-Neck2)")
    check(re.search(r'PBRADII="\$\{PBRADII:-4\}"', txt),
          "PBRadii 默认应为 4 (mbondi3), 与 igb=8 配套")


# =========================================================================
# 4. JSON 阶段计划
# =========================================================================
@case("plan: assets/plan/*.json 结构合法且引用的 mdp 存在")
def t_plans():
    pdir = SKILL / "assets/plan"
    check(pdir.is_dir(), "缺少 assets/plan 目录")
    files = sorted(pdir.glob("*.json"))
    check(files, "assets/plan 里没有 json 计划")
    for f in files:
        obj = json.loads(read(f))
        stages = obj.get("stages") if isinstance(obj, dict) else obj
        check(stages, f"{f.name}: 没有 stages")
        for st in stages:
            for key in ("name", "mdp"):
                check(key in st, f"{f.name}: stage 缺 {key} 字段")
            check(isinstance(st["name"], str) and st["name"],
                  f"{f.name}: stage 的 name 必须是非空字符串")
            mdp = st["mdp"]
            cands = [SKILL / mdp, pdir / mdp, SKILL / "assets/mdp" / mdp]
            check(any(c.is_file() for c in cands),
                  f"{f.name}: stage {st['name']} 引用的 mdp 不存在: {mdp}")


@case("plan: 阶段名唯一(否则会互相覆盖产物)")
def t_plan_unique_names():
    for f in sorted((SKILL / "assets/plan").glob("*.json")):
        obj = json.loads(read(f))
        stages = obj.get("stages") if isinstance(obj, dict) else obj
        names = [s["name"] for s in stages]
        check(len(names) == len(set(names)),
              f"{f.name}: 阶段名重复 {names}")


# =========================================================================
# 5. skill 结构与文档一致性
# =========================================================================
@case("skill: SKILL.md 的 frontmatter name 与目录名一致")
def t_frontmatter():
    txt = read(SKILL / "SKILL.md")
    m = re.search(r"^name:\s*(\S+)\s*$", txt, re.M)
    check(m, "SKILL.md 缺少 name 字段")
    check(m.group(1) == SKILL.name,
          f"SKILL.md 的 name({m.group(1)}) 与目录名({SKILL.name}) 不一致")
    check(re.search(r"^description:\s*\S", txt, re.M), "SKILL.md 缺少 description")


@case("skill: SKILL.md 参考文件列表里的路径都存在")
def t_refs_exist():
    txt = read(SKILL / "SKILL.md")
    paths = set(re.findall(r"`((?:scripts|references|assets)/[^`]+)`", txt))
    missing = []
    for p in sorted(paths):
        if "*" in p:
            if not list(SKILL.glob(p)):
                missing.append(p)
        elif not (SKILL / p).exists():
            missing.append(p)
    check(not missing, "SKILL.md 引用了不存在的文件: " + ", ".join(missing))


@case("skill: README 里提到的脚本都存在")
def t_readme_scripts():
    txt = read(SKILL.parent.parent / "README.md") if (SKILL.parent.parent / "README.md").is_file() else ""
    if not txt:
        return  # 仓库根的 README 不在 skill 目录里, 跳过
    paths = set(re.findall(r"`?scripts/([A-Za-z0-9_.-]+\.(?:sh|py))`?", txt))
    missing = [p for p in sorted(paths) if not (SKILL / "scripts" / p).is_file()]
    check(not missing, "README 引用了不存在的脚本: " + ", ".join(missing))


def _repo_root():
    """定位仓库根(含 LICENSE 的那一层)。"""
    env = os.environ.get("SKILL_REPO_ROOT")
    if env and Path(env).is_dir():
        return Path(env)
    p = SKILL.parent.parent
    if (p / "LICENSE").is_file():
        return p
    return None


@case("skill: 仓库根有 LICENSE 与 NOTICE")
def t_license():
    root = _repo_root()
    if root is None:
        print("      (跳过: 找不到仓库根, 设 SKILL_REPO_ROOT 指向仓库即可检查)")
        return
    check((root / "LICENSE").is_file(), f"{root} 缺少 LICENSE")
    txt = read(root / "LICENSE")
    check("GNU GENERAL PUBLIC LICENSE" in txt.upper(), "LICENSE 不是 GPL 文本")
    check((root / "NOTICE").is_file(),
          f"{root} 缺少 NOTICE(第三方组件与衍生关系说明)")


@case("skill: references 里的 7 篇文档齐全")
def t_references():
    need = ["analysis.md", "docking.md", "equilibration.md", "force-fields.md",
            "gromacs-2025-api.md", "mmpbsa.md", "system-setup.md"]
    have = sorted(p.name for p in (SKILL / "references").glob("*.md"))
    for n in need:
        check(n in have, f"缺少 references/{n}")


# =========================================================================
# 6. 隐私与许可
# =========================================================================
@case("隐私: 文件里不含本机用户路径 / token / 私钥")
def t_no_secrets():
    pats = {
        "token": r"gh[pousr]_[A-Za-z0-9]{16,}",
        "私钥": r"-----BEGIN [A-Z ]*PRIVATE KEY-----",
        "API key": r"(?i)(api[_-]?key|secret[_-]?key)\s*[=:]\s*['\"][A-Za-z0-9]{16,}",
    }
    hits = []
    for f in SKILL.rglob("*"):
        if not f.is_file() or f.suffix in (".png", ".pdf"):
            continue
        try:
            txt = read(f)
        except Exception:
            continue
        for label, pat in pats.items():
            if re.search(pat, txt):
                hits.append(f"{f.name}: {label}")
    check(not hits, "发现疑似敏感信息: " + ", ".join(hits))


@case("隐私: 不含硬编码的本机用户目录")
def t_no_local_paths():
    hits = []
    for f in SKILL.rglob("*"):
        if not f.is_file():
            continue
        try:
            txt = read(f)
        except Exception:
            continue
        for m in re.finditer(r"/home/[A-Za-z0-9_.-]+/", txt):
            hits.append(f"{f.name}: {m.group(0)}")
    check(not hits, "发现硬编码的用户目录: " + ", ".join(sorted(set(hits))))


@case("许可: 不集成 AGPL-3.0 的 charmm2gmx 脚本")
def t_no_agpl_code():
    hits = []
    for f in (SKILL / "scripts").glob("*"):
        if not f.is_file():
            continue
        txt = read(f)
        if "charmm2gmx" in txt and ("GNU Affero" in txt or "AGPL" in txt):
            hits.append(f.name)
    check(not hits, "scripts/ 里不应包含 AGPL-3.0 的代码: " + ", ".join(hits))


# =========================================================================
# 7. 需要真跑 GROMACS 的慢测试
# =========================================================================
@case("实跑: mdp 模板 grompp 零警告(需要现成体系)")
def t_grompp_clean(fast=False):
    work = Path.home() / "mdflow_test"
    if not (work / "solv_ions.gro").is_file():
        print("      (跳过: 没有 ~/mdflow_test 体系, 用 mdflow.py prepare+solvate 造一个)")
        return
    for m in ("em", "nvt", "npt", "md"):
        f = SKILL / "assets/mdp" / f"{m}.mdp"
        p = subprocess.run(
            ["gmx", "grompp", "-f", str(f), "-c", str(work / "solv_ions.gro"),
             "-p", str(work / "complex.top"), "-o", "/tmp/_rt.tpr"],
            capture_output=True, text=True, cwd=str(work))
        out = p.stdout + p.stder
        nwarn = out.count("\nWARNING") + out.count("\nERROR")
        check(p.returncode == 0 and nwarn == 0,
              f"{m}.mdp: exit={p.returncode} 警告/错误={nwarn}")


# =========================================================================
def main():
    ap = argparse.ArgumentParser(description="skill 回归测试")
    ap.add_argument("--fast", action="store_true", help="跳过需要跑 GROMACS 的慢用例")
    ap.add_argument("-k", default=None, help="只跑名字里含该关键字的用例")
    a = ap.parse_args()

    cases = []
    for name, obj in sorted(globals().items()):
        if callable(obj) and hasattr(obj, "_case_name"):
            if a.k and a.k not in obj._case_name:
                continue
            cases.append(obj)

    print(f"{BOLD}protein-ligand-md-workflow 回归测试{NC}  ({len(cases)} 个用例)")
    print(f"  skill: {SKILL}\n")
    failed = 0
    for fn in cases:
        name = fn._case_name
        try:
            if fn.__code__.co_argcount:
                fn(fast=a.fast)
            else:
                fn()
            print(f"  {GREEN}[PASS]{NC} {name}")
        except AssertionError as exc:
            failed += 1
            print(f"  {RED}[FAIL]{NC} {name}")
            print(f"         {exc}")
        except Exception as exc:
            failed += 1
            print(f"  {RED}[ERR ]{NC} {name}")
            print(f"         {type(exc).__name__}: {exc}")
    print()
    total = len(cases)
    print(f"{BOLD}结果:{NC} {total - failed}/{total} 通过"
          + (f", {failed} 失败" if failed else ", 全部通过"))
    return 1 if failed else 0



# =========================================================================
# 8. 报告与作图的硬性约束(这两处出过真实故障, 加回归守卫)
# =========================================================================
@case("报告: make_report.py 具备渲染界面分析页所需的全部要素")
def t_report_interface_page():
    txt = read(SKILL / "scripts" / "make_report.py")
    check(re.search(r"^import json\s*$", txt, re.M),
          "make_report.py 缺少 import json(会让界面分析静默丢失)")
    check("st_h2 = ParagraphStyle(" in txt,
          "make_report.py 未定义 st_h2 样式(渲染界面残基表会崩溃)")
    check(re.search(r'\("interface_map\.png",\s*"interface_map\.png"', txt),
          "FIGURES 里 interface_map.png 的 png 字段是 None, 界面热图不会进报告")
    check(re.search(r'\("dssp\.png",\s*"dssp\.png"', txt),
          "FIGURES 里 dssp.png 的 png 字段是 None, 二级结构图不会进报告")


@case("分析: Ramachandran 用散点图, 脚本路径不依赖调用者")
def t_rama_scatter():
    txt = read(SKILL / "scripts" / "analyze.sh")
    check("SCRIPT_DIR=" in txt, "analyze.sh 没有定义 SCRIPT_DIR")
    check("plot_rama.py" in txt, "analyze.sh 没有调用 plot_rama.py")
    check(re.search(r'python3\s+"?\$SCRIPT_DIR/plot_rama\.py"?', txt),
          "analyze.sh 里 plot_rama.py 的路径没用 $SCRIPT_DIR, 调用时会找不到脚本")
    check(not re.search(r'python3\s+"/scripts/plot_rama\.py"', txt),
          "analyze.sh 里 plot_rama.py 的路径写成了 /scripts/... 绝对路径")



# =========================================================================
# 9. 报告内容完整性与排版
# =========================================================================
def _load_make_report():
    import importlib.util
    path = SKILL / "scripts" / "make_report.py"
    spec = importlib.util.spec_from_file_location("_mr_under_test", str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@case("报告: MM/PBSA 的 ΔG 能从 gmx_MMPBSA 输出里正确解析")
def t_mmpbsa_parse():
    import tempfile
    mod = _load_make_report()
    # 1.7 格式: Average | SD(Prop.) | SD | SEM(Prop.) | SEM | Block SD | Block SEM
    new_fmt = (
        "Delta (Complex - Receptor - Ligand):\n"
        "Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)"
        "        SEM   Block SD  Block SEM\n"
        + "-" * 100 + "\n"
        "\u0394TOTAL                   17.25    1740390.66      50.49"
        "    524747.53      15.22       2.35       1.66\n"
    )
    # 老格式: Value | Average | Std Dev | Std Err of Est. | Number of Samples
    old_fmt = (
        "Delta (Complex - Receptor - Ligand):\n"
        "Energy Component       Value      Average      Std Dev"
        "    Std Err of Est.   Number of Samples\n"
        + "-" * 100 + "\n"
        "\u0394TOTAL                 -17.25     -17.25       1.70"
        "       1.66              11\n"
    )
    for text, want_dg, want_bs in ((new_fmt, 17.25, 1.66),
                                   (old_fmt, -17.25, 1.66)):
        with tempfile.TemporaryDirectory() as td:
            Path(td, "mmpbsa.out").write_text(text, encoding="utf-8")
            got = mod._parse_mmpbsa(Path(td))
        check(got is not None, "解析 mmpbsa.out 失败")
        check(abs(got["dG"] - want_dg) < 1e-6,
              "ΔG 解析错误: 期望 %s, 实际 %s" % (want_dg, got["dG"]))
        check(got.get("block_sem") is not None
              and abs(got["block_sem"] - want_bs) < 1e-6,
              "Block SEM 解析错误: 期望 %s, 实际 %s" % (want_bs, got.get("block_sem")))


@case("报告: 正文字号为 12 号")
def t_body_font_size():
    txt = read(SKILL / "scripts" / "make_report.py")
    m = re.search(r"^BODY_FS\s*=\s*([\d.]+)", txt, re.M)
    check(m, "make_report.py 里找不到 BODY_FS 常量")
    check(float(m.group(1)) == 12.0,
          "正文字号应为 12 号, 实际 %s" % m.group(1))
    check("BODY_FS" in txt and "fontSize=BODY_FS" in txt,
          "正文样式没有走 BODY_FS 常量")


@case("报告: 界面分析/能量图/产物清单都在报告里")
def t_report_sections():
    txt = read(SKILL / "scripts" / "make_report.py")
    for key, why in (("\u754c\u9762\u76f8\u4e92\u4f5c\u7528\u5206\u6790", "界面分析页"),
                     ("RIN", "残基-残基接触网络表"),
                     ("\u5206\u6790\u4ea7\u7269\u6e05\u5355", "分析产物清单页"),
                     ("energy.png", "体系能量分量图")):
        check(key in txt, "报告里缺少" + why)



if __name__ == "__main__":
    sys.exit(main())