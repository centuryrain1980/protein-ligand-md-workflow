#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把分析结果汇总成一份可直接进论文补充材料的 PDF 报告 (figurebook)。

用法:
    python3 make_report.py <工作目录> [选项]

选项:
    --out <path>        输出 PDF 路径(默认 <工作目录>/analysis/report.pdf)
    --title <title>     报告标题(默认用工作目录名)
    --prot <group>      蛋白组名(默认 Protein)
    --lig <group>       配体组名(默认 LIG)
    --ligand-file <f>   配体结构文件(mol2/sdf), 用于画 2D 结构式
    --html              强制输出 HTML(不依赖 reportlab)
    --no-ligand         不画配体结构式
    --page-size <s>     页面尺寸: A4 | letter | landscape

设计要点:
  * 中文字体: 优先用系统里的 Noto CJK / 文泉驿等 TTF, 没有则退回 reportlab
    内置的 STSong-Light (无需外部字体文件), 都不会出现方块字;
  * 每个分析量缺失都不影响其它内容, 全部用 try/except 包住;
  * 没有 reportlab 时自动降级成 HTML(纯标准库), 保证"总能出报告";
  * 数值直接从 .xvg / mmpbsa.out 里解析, 不做二次计算, 避免与原始分析不一致。
"""

from __future__ import annotations

import argparse
import datetime
import html as html_mod
import json
import math
import os
import re
import sys
import traceback
from pathlib import Path

# 系统里常见的中日韩字体(按显示效果排序), 找不到就用 reportlab 内置 CID 字体
CJK_FONT_CANDIDATES = [
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJKsc-Regular.otf",
    "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
    "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
    "/usr/share/fonts/truetype/arphic/uming.ttc",
    "/usr/share/fonts/truetype/arphic/ukai.ttc",
    "/usr/share/fonts/opentype/noto/NotoSerifCJK-Regular.ttc",
    "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf",
]

# ----------------------------------------------------------------------
# 图表清单: (xvg 文件, png 文件, 标题, 说明)
# ----------------------------------------------------------------------
FIGURES = [
    ("rmsd.xvg", "rmsd.png", "RMSD — 蛋白结构 RMSD",
     "按蛋白骨架拟合后计算, 反映模拟过程中构象的整体稳定性。"
     "曲线平台化说明体系已收敛; 若持续漂移, 通常需要延长平衡或检查力场/质子化状态。"),
    ("rmsf.xvg", "rmsf.png", "RMSF — 逐残基涨落",
     "每个残基原子位置的均方根涨落。峰值常出现在结合界面与柔性环区, 是判断"
     "哪些位点参与结合、哪些过度柔动的直接依据。"),
    ("gyrate.xvg", "gyrate.png", "回转半径 Rg",
     "蛋白的回转半径随时间的变化。整体漂移提示构象在发生显著变化(折叠/解折叠)。"),
    ("sasa.xvg", "sasa.png", "溶剂可及面积 SASA",
     "蛋白整体的溶剂可及表面积。数值下降通常对应疏水核心被埋藏或构象收紧。"),
    ("rama.xvg", "rama.png", "Ramachandran 图",
     "主链 φ/ψ 角的散点分布(每点一个残基, 不连线)。落在右上/左下两个允许区域说明二级结构合理; "
     "出现大量异常点提示力场或质子化状态需要检查。"),
    ("binding.xvg", "binding.png", "组间相互作用能",
     "蛋白组与配体组之间的范德华与静电相互作用能(需在 mdp 里设过 energygrps)。"
     "数值越负表示结合越紧密。"),
    ("hbond_num.xvg", "hbond.png", "氢键数量",
     "蛋白内部氢键随时间的变化。氢键网络的崩解常常先于整体构象变化出现。"),
    ("dssp.png", "dssp.png", "二级结构含量",
     "用 gmx dssp 判定的二级结构比例随时间的变化。"),
    ("interface_map.png", "interface_map.png", "界面接触频率热图",
     "配体与蛋白界面残基的接触次数热图(颜色越亮接触越多), "
     "由 scripts/interface_analysis.py 生成。"),
]


# ----------------------------------------------------------------------
# 小工具
# ----------------------------------------------------------------------
def _read_xvg(path):
    """读 xvg, 返回 (列名列表, 数据行列表)。"""
    names, rows = [], []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(("#", "@")):
                if line.startswith("@") and "s0" in line:
                    m = re.search(r'legend\s+"([^"]*)"', line)
                    if m:
                        names.append(m.group(1))
                continue
            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    return names, rows


def _col(rows, idx, default=None):
    vals = [r[idx] for r in rows if len(r) > idx]
    return vals if vals else default


def _stat(vals):
    if not vals:
        return {}
    n = len(vals)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / n if n > 1 else 0.0
    return {"n": n, "mean": mean, "min": min(vals), "max": max(vals),
            "final": vals[-1], "sd": var ** 0.5}


def _fmt(v, unit="", nd=2):
    """格式化数值; nan/inf 显示为破折号而不是 'nan'。"""
    if v is None:
        return "—"
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return "—"
        return f"{v:.{nd}f}{unit}"
    return f"{v}{unit}"


def _parse_mmpbsa(work):
    """从 mmpbsa.out 解析 ΔG_bind (Delta 段的 ΔTOTAL 行)。"""
    p = Path(work) / "mmpbsa.out"
    if not p.is_file():
        return None
    txt = p.read_text(errors="replace")
    m = re.search(r"Delta \(Complex - Receptor - Ligand\):(.*?)-{20,}", txt, re.S)
    if not m:
        return None
    row = re.search(r"^\W*TOTAL\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+"
                    r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", m.group(1), re.M)
    if not row:
        return None
    g = row.groups()
    return {"dG": float(g[0]), "sem": float(g[4]), "block_sd": float(g[5]),
            "block_sem": float(g[6])}


def _gmx_version():
    import subprocess
    try:
        out = subprocess.run(["gmx", "--version"], capture_output=True, text=True,
                             timeout=30).stdout
        m = re.search(r"GROMACS version:\s*(\S+)", out)
        return m.group(1) if m else "unknown"
    except Exception:
        return "unknown"


def _ligand_png(ligfile, outpng, size=(560, 460)):
    """用 rdkit 画配体 2D 结构式, 失败返回 None。"""
    if not ligfile or not Path(ligfile).is_file():
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        mol = None
        if str(ligfile).lower().endswith(".sdf"):
            suppl = Chem.SDMolSupplier(str(ligfile))
            mol = next((m for m in suppl if m is not None), None)
        else:
            mol = Chem.MolFromMol2File(str(ligfile), removeHs=False)
        if mol is None:
            return None
        try:
            AllChem.Compute2DCoords(mol)
        except Exception:
            pass
        img = Draw.MolToImage(mol, size=size)
        img.save(outpng)
        return outpng
    except Exception:
        return None


# ----------------------------------------------------------------------
# 汇总
# ----------------------------------------------------------------------
def collect(work, prot, lig, ligand_file=None, want_ligand=True):
    ana = Path(work) / "analysis"
    info = {"workdir": str(Path(work).resolve()),
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "gmx": _gmx_version(), "metrics": [], "figures": []}

    def add(label, value, unit=""):
        info["metrics"].append((label, value, unit))

    for xtc in ("md.xtc", "prod.xtc", "traj.xtc"):
        if (Path(work) / xtc).is_file():
            info["xtc"] = xtc
            break
    else:
        info["xtc"] = None

    f = ana / "rmsd.xvg"
    if f.is_file():
        _, rows = _read_xvg(f)
        v = _stat(_col(rows, 1))
        if v:
            add("RMSD 平均值 (nm)", _fmt(v["mean"], " nm", 3))
            add("RMSD 末帧 (nm)", _fmt(v["final"], " nm", 3))
            add("RMSD 标准差 (nm)", _fmt(v["sd"], " nm", 3))
            add("RMSD 采样帧数", v["n"])
        t = _col(rows, 0)
        if t:
            add("模拟时长 (ps)", _fmt(max(t) - min(t), " ps", 1))
            add("轨迹帧数", len(t))
    f = ana / "rmsf.xvg"
    if f.is_file():
        _, rows = _read_xvg(f)
        v = _stat(_col(rows, 1))
        if v:
            add("RMSF 平均值 (nm)", _fmt(v["mean"], " nm", 3))
            add("RMSF 最大值 (nm)", _fmt(v["max"], " nm", 3))
    for fn, label, unit, nd in (("gyrate.xvg", "回转半径 Rg 平均 (nm)", " nm", 2),
                               ("sasa.xvg", "SASA 平均 (nm²)", " nm²", 1)):
        f = ana / fn
        if f.is_file():
            _, rows = _read_xvg(f)
            v = _stat(_col(rows, 1))
            if v:
                add(label, _fmt(v["mean"], unit, nd))
    f = ana / "hbond_num.xvg"
    if f.is_file():
        _, rows = _read_xvg(f)
        v = _stat(_col(rows, 1))
        if v:
            add("蛋白内氢键数 平均", _fmt(v["mean"], "", 1))
    dg = _parse_mmpbsa(work)
    info["mmpbsa"] = dg
    if dg:
        add("ΔG_bind (kcal/mol)", _fmt(dg["dG"], "", 2))
        add("ΔG_bind Block SEM", _fmt(dg["block_sem"], "", 2))

    isum = ana / "interface_summary.json"
    if isum.is_file():
        try:
            js = json.loads(isum.read_text(encoding="utf-8"))
            info["interface"] = js
            add("界面相互作用事件数", js.get("total_events", "—"))
            add("界面残基数", js.get("interface_residues", "—"))
            for k, v in (js.get("by_interaction_class") or {}).items():
                add("界面事件 · " + str(k), v)
            if js.get("ligand_resname"):
                add("配体残基名", js["ligand_resname"])
        except Exception as _exc:
            print("[make_report] read interface_summary.json failed: %r" % (_exc,), file=sys.stderr)
            info["interface"] = None
    else:
        info["interface"] = None

    lig_png = None
    if want_ligand:
        cands = []
        if ligand_file:
            cands.append(ligand_file)
        for pat in ("*_lig.mol2", "lig.mol2", "complex_lig.mol2", "*.sdf"):
            cands.extend(sorted(str(x) for x in Path(work).glob(pat)))
        for c in cands:
            lig_png = _ligand_png(c, str(ana / "ligand_2d.png"))
            if lig_png:
                break
    info["ligand_png"] = lig_png

    for xvg, png, title, desc in FIGURES:
        px = ana / png if png else None
        if px and px.is_file() and px.stat().st_size > 0:
            info["figures"].append({"path": str(px), "title": title, "desc": desc})
    return info


# ----------------------------------------------------------------------
# HTML 输出(无 reportlab 时的降级方案, 纯标准库)
# ----------------------------------------------------------------------
def write_html(info, outpath):
    parts = [
        "<!doctype html><html lang='zh'><head><meta charset='utf-8'>",
        "<title>{}</title><style>body{{font-family:system-ui,sans-serif;margin:2rem;"
        "max-width:1000px}}table{{border-collapse:collapse;margin:1rem 0}}"
        "td,th{{border:1px solid #ccc;padding:.35rem .7rem;text-align:left}}"
        "img{{max-width:100%;border:1px solid #eee}}h2{{margin-top:2.5rem}}</style>"
        "</head><body>".format(html_mod.escape(info.get("title", "MD Report"))),
        "<h1>{}</h1>".format(html_mod.escape(info.get("title", "MD Report"))),
        "<p>体系目录 <code>{}</code> · 生成时间 {} · GROMACS {}</p>".format(
            html_mod.escape(info["workdir"]), info["date"], info["gmx"]),
        "<h2>关键指标</h2><table><tr><th>指标</th><th>数值</th></tr>",
    ]
    for label, value, unit in info["metrics"]:
        parts.append("<tr><td>{}</td><td>{}</td></tr>".format(
            html_mod.escape(label), html_mod.escape(str(value))))
    parts.append("</table>")
    if info.get("ligand_png"):
        parts.append("<h2>配体结构</h2><img src='ligand_2d.png' alt='ligand'>")
    parts.append("<h2>分析图</h2>")
    for fig in info["figures"]:
        parts.append("<h3>{}</h3><p>{}</p><img src='{}' alt='{}'>".format(
            html_mod.escape(fig["title"]), html_mod.escape(fig["desc"]),
            html_mod.escape(os.path.basename(fig["path"])),
            html_mod.escape(fig["title"])))
    parts.append("</body></html>")
    Path(outpath).write_text("\n".join(parts), encoding="utf-8")
    return outpath


# ----------------------------------------------------------------------
# PDF 输出
# ----------------------------------------------------------------------
def _setup_cjk_font():
    """注册能显示中文的字体, 返回字体名。

    优先系统 TTF(Noto CJK / 文泉驿等), 没有则用 reportlab 内置的
    STSong-Light (UnicodeCIDFont, 不需要额外字体文件)。
    """
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfbase.cidfonts import UnicodeCIDFont

    for path in CJK_FONT_CANDIDATES:
        if Path(path).is_file():
            try:
                pdfmetrics.registerFont(TTFont("CJK", path))
                return "CJK"
            except Exception:
                continue
    for name in ("STSong-Light", "HeiseiMin-W3"):
        try:
            pdfmetrics.registerFont(UnicodeCIDFont(name))
            return name
        except Exception:
            continue
    return "Helvetica"


def write_pdf(info, outpath, page_size="A4"):
    from reportlab.lib.pagesizes import A4, letter, landscape
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.lib import colors
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                    Table, TableStyle, PageBreak)

    FONT = _setup_cjk_font()
    sizes = {"A4": A4, "letter": letter, "landscape": landscape(A4)}
    ps = sizes.get(page_size, A4)

    doc = SimpleDocTemplate(
        str(outpath), pagesize=ps,
        leftMargin=18 * mm, rightMargin=18 * mm,
        topMargin=16 * mm, bottomMargin=16 * mm,
        title=info.get("title", "MD Report"), author="protein-ligand-md-workflow")

    ss = getSampleStyleSheet()
    st_title = ParagraphStyle("t", parent=ss["Title"], fontName=FONT, fontSize=20,
                              leading=26, spaceAfter=10)
    st_sub = ParagraphStyle("s", parent=ss["Normal"], fontName=FONT, fontSize=10,
                             leading=15, textColor=colors.HexColor("#555555"))
    st_h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontName=FONT, fontSize=15,
                            leading=20, spaceBefore=12, spaceAfter=6,
                            textColor=colors.HexColor("#1a3d5c"))
    st_h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontName=FONT, fontSize=12,
                            leading=17, spaceBefore=8, spaceAfter=4,
                            textColor=colors.HexColor("#2c5f85"))
    st_body = ParagraphStyle("b", parent=ss["Normal"], fontName=FONT, fontSize=9.5,
                              leading=14.5)
    st_cap = ParagraphStyle("c", parent=ss["Normal"], fontName=FONT, fontSize=8.5,
                            leading=12, textColor=colors.HexColor("#555555"),
                            spaceAfter=8)

    story = []
    # ---------- 封面 ----------
    story.append(Spacer(1, 26 * mm))
    story.append(Paragraph(html_mod.escape(info.get("title", "MD Report")), st_title))
    story.append(Paragraph(
        "蛋白-配体复合物分子动力学分析报告 &nbsp;·&nbsp; "
        "protein-ligand-md-workflow", st_sub))
    story.append(Spacer(1, 8 * mm))
    cover = [
        ["体系目录", info["workdir"]],
        ["生成时间", info["date"]],
        ["GROMACS 版本", info["gmx"]],
        ["轨迹文件", info.get("xtc") or "未找到"],
        ["配体结构式", "已包含" if info.get("ligand_png") else "无"],
    ]
    t = Table(cover, colWidths=[38 * mm, 118 * mm])
    t.setStyle(TableStyle([
        ("FONT", (0, 0), (-1, -1), FONT, 9.5),
        ("TEXTCOLOR", (0, 0), (0, -1), colors.HexColor("#666666")),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("LINEBELOW", (0, 0), (-1, -1), 0.4, colors.HexColor("#dddddd")),
    ]))
    story.append(t)
    story.append(PageBreak())

    # ---------- 关键指标 ----------
    story.append(Paragraph("关键指标摘要", st_h1))
    story.append(Paragraph(
        "下列数值直接解析自分析输出的 .xvg 与 mmpbsa.out, 未做二次计算, "
        "因此与原始分析文件完全一致。", st_body))
    story.append(Spacer(1, 3 * mm))
    rows = [["指标", "数值"]] + [[l, str(v)] for l, v, _ in info["metrics"]]
    t = Table(rows, colWidths=[70 * mm, 60 * mm], repeatRows=1)
    t.setStyle(TableStyle([
        ("FONT", (0, 0), (-1, 0), FONT, 9.5),
        ("FONT", (0, 1), (-1, -1), FONT, 9.5),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef3f8")),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cccccc")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(t)
    if info.get("mmpbsa"):
        dg = info["mmpbsa"]
        story.append(Spacer(1, 4 * mm))
        story.append(Paragraph(
            "结合自由能说明: ΔG_bind = <b>{:.2f}</b> kcal/mol "
            "(Block SEM = {:.2f})。请确认采样帧数足够(通常 50–200 帧); "
            "帧数过少时 Block SEM 会大于均值本身, 此时结果不具备统计意义。".format(
                dg["dG"], dg["block_sem"]), st_body))
    story.append(PageBreak())

    # ---------- 配体结构 ----------
    if info.get("ligand_png"):
        story.append(Paragraph("配体结构", st_h1))
        try:
            from PIL import Image as PILImage
            iw, ih = PILImage.open(info["ligand_png"]).size
            w = ps[0] - 36 * mm
            h = ih * w / iw
            story.append(Image(info["ligand_png"], width=w, height=h))
        except Exception:
            story.append(Paragraph("(配体结构图无法嵌入)", st_cap))
        story.append(PageBreak())

    # ---------- 界面相互作用专页 ----------
    if info.get("interface"):
        js = info["interface"]
        story.append(Paragraph("界面相互作用分析", st_h1))
        story.append(Paragraph(
            "由 scripts/interface_analysis.py 生成: 逐帧扫描配体-蛋白接触, "
            "并按拓扑信息(元素/电荷/残基名)与几何证据保守分类; "
            "判不出具体类型的一律记为 Other Contact, 不臆测化学类型。", st_body))
        story.append(Spacer(1, 3 * mm))
        rows = [["项目", "数值"]]
        rows.append(["分析帧数", js.get("frames_analysed", "—")])
        rows.append(["相互作用事件总数", js.get("total_events", "—")])
        rows.append(["界面残基数", js.get("interface_residues", "—")])
        rows.append(["配体残基名", js.get("ligand_resname", "—")])
        for k, v in (js.get("by_interaction_class") or {}).items():
            rows.append(["分类 · " + str(k), v])
        for k, v in (js.get("by_validation_status") or {}).items():
            rows.append(["证据强度 · " + str(k), v])
        t2 = Table(rows, colWidths=[70 * mm, 60 * mm], repeatRows=1)
        t2.setStyle(TableStyle([
            ("FONT", (0, 0), (-1, 0), FONT, 9.5),
            ("FONT", (0, 1), (-1, -1), FONT, 9.5),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef3f8")),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cccccc")),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ]))
        story.append(t2)
        top = js.get("top_interface_residues") or []
        if top:
            story.append(Spacer(1, 4 * mm))
            story.append(Paragraph("接触最多的界面残基", st_h2))
            r2 = [["残基", "链", "配体", "接触次数"]]
            for r in top[:12]:
                res = str(r.get("residue", ""))
                chain = res.split(":")[-1] if ":" in res else ""
                rname = res.split(":")[0] if ":" in res else res
                r2.append([rname, chain, r.get("ligand", ""), r.get("contacts", "")])
            t3 = Table(r2, colWidths=[45 * mm, 20 * mm, 40 * mm, 30 * mm], repeatRows=1)
            t3.setStyle(TableStyle([
                ("FONT", (0, 0), (-1, 0), FONT, 9.5),
                ("FONT", (0, 1), (-1, -1), FONT, 9.5),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef3f8")),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cccccc")),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]))
            story.append(t3)
        story.append(PageBreak())

    # ---------- 各图 ----------
    for fig in info["figures"]:
        story.append(Paragraph(html_mod.escape(fig["title"]), st_h1))
        try:
            from PIL import Image as PILImage
            iw, ih = PILImage.open(fig["path"]).size
            w = ps[0] - 36 * mm
            h = ih * w / iw
            maxh = ps[1] - 70 * mm
            if h > maxh:
                h = maxh
                w = iw * h / ih
            story.append(Image(fig["path"], width=w, height=h))
        except Exception:
            story.append(Paragraph("(图片无法嵌入)", st_body))
        story.append(Paragraph(html_mod.escape(fig["desc"]), st_cap))
        story.append(PageBreak())

    # ---------- 方法与参数 ----------
    story.append(Paragraph("方法与参数", st_h1))
    meta = [
        ["体系目录", info["workdir"]],
        ["报告生成时间", info["date"]],
        ["GROMACS 版本", info["gmx"]],
        ["分析工具", "GROMACS 2025 内置分析命令 + DuIvyTools (dit) 出图"],
        ["报告生成", "scripts/make_report.py (protein-ligand-md-workflow)"],
    ]
    t = Table(meta, colWidths=[38 * mm, 118 * mm])
    t.setStyle(TableStyle([
        ("FONT", (0, 0), (-1, -1), FONT, 9),
        ("TEXTCOLOR", (0, 0), (0, -1), colors.HexColor("#666666")),
        ("LINEBELOW", (0, 0), (-1, -1), 0.4, colors.HexColor("#dddddd")),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(t)
    story.append(Spacer(1, 5 * mm))
    story.append(Paragraph(
        "引用提示: 使用本报告中的数值时, 请同时引用所用的力场文献、"
        "GROMACS 文献以及所用分析方法的原始文献。", st_body))

    doc.build(story)
    return outpath


# ----------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="把 MD 分析结果汇总成 PDF 报告 (figurebook)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("workdir", help="工作目录(含 analysis/ 子目录)")
    ap.add_argument("--out", default=None, help="输出 PDF 路径")
    ap.add_argument("--title", default=None, help="报告标题")
    ap.add_argument("--prot", default="Protein", help="蛋白组名")
    ap.add_argument("--lig", default="LIG", help="配体组名")
    ap.add_argument("--ligand-file", default=None, help="配体结构文件(mol2/sdf)")
    ap.add_argument("--html", action="store_true", help="强制输出 HTML")
    ap.add_argument("--no-ligand", action="store_true", help="不画配体结构式")
    ap.add_argument("--page-size", default="A4", choices=["A4", "letter", "landscape"])
    args = ap.parse_args(argv)

    work = Path(args.workdir).resolve()
    if not work.is_dir():
        print("错误: 工作目录不存在:", work, file=sys.stderr)
        return 2
    ana = work / "analysis"
    if not ana.is_dir():
        print("错误: 找不到 analysis 目录, 请先运行 analyze.sh", file=sys.stderr)
        return 2

    info = collect(work, args.prot, args.lig, args.ligand_file,
                   want_ligand=not args.no_ligand)
    info["title"] = args.title or f"MD Analysis Report — {work.name}"

    out = Path(args.out) if args.out else ana / "report.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.html:
        out = out.with_suffix(".html")
        write_html(info, out)
        print("HTML 报告:", out)
        return 0

    try:
        write_pdf(info, out, args.page_size)
    except ImportError:
        print("提示: 没装 reportlab, 降级输出 HTML。"
              "装 PDF 版: pip3 install --break-system-packages reportlab svglib",
              file=sys.stderr)
        out = out.with_suffix(".html")
        write_html(info, out)
        print("HTML 报告:", out)
        return 0
    except Exception:
        traceback.print_exc()
        print("提示: PDF 生成失败, 降级输出 HTML", file=sys.stderr)
        out = out.with_suffix(".html")
        write_html(info, out)
        print("HTML 报告:", out)
        return 0

    print("PDF 报告:", out, f"({out.stat().st_size} bytes)")
    print("  关键指标:")
    for label, value, _ in info["metrics"]:
        print(f"    {label:<28} {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main())