# protein-ligand-md-workflow

> **蛋白-配体复合物的分子动力学全流程 Codex Skill**
> 从 PDB 结构和配体 SMILES 出发，一条流水线走完：**分子对接 → 体系建模 → 力场选型 → 参数化 → 溶剂化 → 平衡 → 生产 MD → MM/PBSA 结合自由能 → 结果分析**。
>
> 底层工具：GROMACS 2025 + DuIvyTools(`dit`) + acpype + gmx_MMPBSA + gnina / AutoDock Vina。
> 所有命令与参数**在 GROMACS 2025 上实测通过**，并把踩过的坑直接写进了配套文档。

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

**English summary** — A Codex skill that automates the end-to-end molecular-dynamics
workflow for protein–ligand complexes: docking, system preparation, force-field
selection, ligand parameterization (GAFF), solvation, equilibration, production MD,
MM/PBSA binding free energies and trajectory analysis. Built and validated against
GROMACS 2025, DuIvyTools, acpype, gmx_MMPBSA 1.7 and gnina 1.3.3.

---

## 目录

- [能做什么](#能做什么)
- [目录结构](#目录结构)
- [环境要求](#环境要求)
- [安装](#安装)
- [快速开始](#快速开始)
- [分模块使用说明](#分模块使用说明)
- [力场选型速查](#力场选型速查)
- [常见问题（踩坑清单）](#常见问题踩坑清单)
- [与 GROMACS 2025 的兼容性说明](#与-gromacs-2025-的兼容性说明)
- [致谢与许可证](#致谢与许可证)

---

## 能做什么

| 阶段 | 能力 | 说明 |
|---|---|---|
| **0. 环境自检** | `setup_env.sh` | 一键检查 GROMACS / dit / acpype / obabel / vina / gnina / 分析库是否就位，并给出安装建议 |
| **0.5 力场选型** | `references/force-fields.md` | 本机可用力场清单、选型决策树、水模型配对表、CHARMM36 官方强制 mdp、AMBER→GROMACS 转换细节 |
| **1. 分子对接** | `dock_ligand.sh` | gnina（Vina 打分 + CNN 深度学习重打分）把配体送进口袋；自动定盒子、转 PDBQT |
| **2. 配体参数化** | `setup_ligand.sh` | **SMILES → 3D → GAFF 参数 → 力场合并 → 复合物组装**，一条命令 |
| **3. 力场合并** | `merge_gaff_ff.py` | 把 GAFF 配体原子类型并入蛋白力场目录（解决 GROMACS 只允许一个 `[ atomtypes ]` 段的问题） |
| **4. 复合物组装** | `assemble_complex.py` | 蛋白 + 配体 → 复合物坐标与拓扑，含残基号 / 原子序号重编 |
| **5. 溶剂化与加离子** | `assets/mdp/ions.mdp` | 建盒、加水（TIP3P/SPC/TIP4P…）、genion 加离子 |
| **6. 平衡** | `assets/mdp/{em,nvt,npt}.mdp` | EM → NVT → NPT，**全部零警告**（GROMACS 2025 把警告当错误） |
| **7. 生产 MD** | `assets/mdp/md.mdp` | 长程轨迹，可选输出组间能量项用于结合能分析 |
| **8. MM/PBSA** | `run_mmpbsa.sh` | gmx_MMPBSA 1.7.0 结合自由能 + 逐残基能量分解，自动汇总结果 |
| **9. 结果分析** | `analyze.sh` | RMSD / RMSF / 能量 / 二级结构 / Ramachandran / 自由能面，批量出图 |

**特色**

- **零警告原则**：所有 mdp 模板按 GROMACS 2025「警告即错误」的规则逐条校正，实测 `grompp` 0 warning。
- **踩坑写进文档**：GROMACS 2025 的 API 变更、`dit ndx_add` 不覆盖已存在文件、gmx_MMPBSA 输出段落误读等，
  全部记录在 `references/` 里，避免重复踩雷。
- **官方文档为依据**：力场与 MM/PBSA 参数全部对照 GROMACS 官方手册与 gmx_MMPBSA 官方文档，并注明出处。
- **可移植**：关键路径均可用环境变量覆盖（`GMXBIN` / `GMXDATA` / `PROT_FF` / `PROT_WATER` / `MMPBSA_HOME`）。

---

## 目录结构

```
protein-ligand-md-workflow/
├── README.md                  ← 本文件
├── LICENSE                    ← GPL-3.0 全文
├── NOTICE                     ← 第三方组件与衍生关系说明
├── .gitignore
└── protein-ligand-md-workflow/    ← skill 本体（安装的就是它）
    ├── SKILL.md               ← skill 主入口：环境、流程总览、最小骨架、硬性约束
    ├── agents/openai.yaml     ← Codex 界面元信息
    ├── assets/mdp/            ← 10 个 mdp 模板（5 个通用 + 5 个 CHARMM36 专用）
    ├── references/            ← 深度参考资料（按需加载，不占上下文）
    │   ├── force-fields.md    ← 力场选型与参数转换
    │   ├── system-setup.md    ← 蛋白/配体准备、组装、溶剂化、离子、索引组
    │   ├── docking.md         ← Vina/gnina 对接、盒子定位、版本对照
    │   ├── equilibration.md   ← EM/NVT/NPT/生产 MD 参数原理与调优
    │   ├── analysis.md        ← 各类分析命令 + dit 命令映射表
    │   ├── mmpbsa.md          ← MM/PBSA 环境、参数、输出格式、踩坑
    │   └── gromacs-2025-api.md← GROMACS 2025 命令参数变更与排错
    └── scripts/               ← 可直接执行的一键脚本
        ├── setup_env.sh       环境自检
        ├── setup_ligand.sh    配体一条命令（SMILES → 复合物）
        ├── setup_mmpbsa.sh    安装 gmx_MMPBSA 环境
        ├── merge_gaff_ff.py   GAFF 力场合并
        ├── assemble_complex.py 复合物组装
        ├── dock_ligand.sh     分子对接
        ├── run_simulation.sh  批量跑平衡 + 生产
        ├── run_mmpbsa.sh      一键 MM/PBSA
        └── analyze.sh         批量分析出图
```

---

## 环境要求

| 组件 | 版本 | 说明 |
|---|---|---|
| **GROMACS** | 2025.x（2025.4 实测） | 必须 2025+；2024 及以前的参数差异很大，脚本不通用 |
| **DuIvyTools** | 0.6.0 | 提供 `dit` 命令，负责 xvg/xpm 可视化与后处理 |
| **acpype** | 2026.x | 小分子配体 GAFF/GAFF2 参数化（需 AmberTools 的 `antechamber`） |
| **AmberTools** | 26.0 | 提供 `antechamber` / `cpptraj` / `tleap` / `sander` |
| **gmx_MMPBSA** | 1.7.0 | MM/PBSA 结合自由能；**要求 Python ≥3.11 且 <3.13** |
| **gnina** | 1.3.3（1.1 亦可） | CNN 重打分对接；直接吃 PDB 受体 |
| **AutoDock Vina** | 1.2.7 | 备选对接引擎 |
| **Open Babel** | 3.1+ | 格式转换（注意 3.x 的 `-o` / `-O` 语义变了） |
| **meeko** | 0.7+ | 受体 PDBQT 制备（走 Vina 路线时需要） |
| **Python 库** | — | `numpy` `scipy` `matplotlib` `pandas` `MDAnalysis` `mdtraj` |

**硬件建议**：8 核以上、16 GB 内存可跑中等体系；生产 MD 建议用 `gmx_mpi` 多进程。
NVIDIA GPU 可显著加速 `mdrun`；对接用 `--device 0` 加速 CNN 打分。

---

## 安装

### 第 1 步：安装 skill

**方式 A —— 用 Codex 自带的 skill 安装器（推荐）**

```bash
python3 ~/.codex/skills/.system/skill-installer/scripts/install-skill-from-github.py \
  --repo centuryrain1980/protein-ligand-md-workflow \
  --path  protein-ligand-md-workflow
```

**方式 B —— 手动克隆**

```bash
git clone https://github.com/centuryrain1980/protein-ligand-md-workflow.git
cp -r protein-ligand-md-workflow/protein-ligand-md-workflow ~/.codex/skills/
```

装好后 **重启 Codex 或新开一轮对话**，skill 才会被加载。

### 第 2 步：安装依赖

在 WSL2 Ubuntu（或 Linux）里：

```bash
# GROMACS 2025（含 AVX2 优化）
# 参考 https://manual.gromacs.org/current/installation.html
# 源码编译时务必带上: -DGMX_SIMD=AVX2_256 -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF

# DuIvyTools
pip3 install --break-system-packages -i https://pypi.tuna.tsinghua.edu.cn/simple DuIvyTools

# 配体参数化与分析库（建议虚拟环境，Ubuntu 26.04 的 rich 包由 Debian 托管）
python3 -m venv ~/venvs/mdtools
source ~/venvs/mdtools/bin/activate
pip install acpype mdanalysis mdtraj meeko rdkit gemmi prody

# Open Babel / Vina / gnina
apt install openbabel
# vina: https://github.com/ccsb-scripps/AutoDock-Vina/releases
# gnina: https://github.com/gnina/gnina/releases
```

### 第 3 步：自检

```bash
bash ~/.codex/skills/protein-ligand-md-workflow/scripts/setup_env.sh
```

会逐项打印 `[OK] / [WARN] / [FAIL]`，缺什么给出对应的安装命令。

### 第 4 步（可选）：MM/PBSA 环境

```bash
bash ~/.codex/skills/protein-ligand-md-workflow/scripts/setup_mmpbsa.sh
```

用 micromamba 建一个独立的 `mmpbsa` 环境（Python 3.12 + AmberTools 26.0 + gmx_MMPBSA 1.7.0），
**不要**装进系统 Python（gmx_MMPBSA 要求 `<3.13`）。

---

## 快速开始

### 最短路径：一条命令从 SMILES 走到复合物

```bash
SKILL=~/.codex/skills/protein-ligand-md-workflow

# 1) 蛋白拓扑
gmx pdb2gmx -f receptor.pdb -o protein.gro -i posre.itp -p topol.top \
            -ff amber99sb -water tip3p -ignh < /dev/null

# 2) 配体：SMILES → 3D → GAFF → 力场合并 → 复合物
bash $SKILL/scripts/setup_ligand.sh ligand.smi topol.top protein.gro complex
```

> `setup_ligand.sh` 的默认力场是 `amber99sb` + `tip3p`，
> 可用 `PROT_FF=oplsaa PROT_WATER=spc` 覆盖。

### 完整流程

```bash
cd ~/projects/<体系名> && mkdir -p mdp && cd mdp
SKILL=~/.codex/skills/protein-ligand-md-workflow

# ── 1. 对接（配体进口袋）────────────────────────────────
bash $SKILL/scripts/dock_ligand.sh receptor.pdb ligand.mol2 30.25 28.91 15.03 22

# ── 2. 蛋白 + 配体（见上面「最短路径」）──────────────────

# ── 3. 溶剂化 ───────────────────────────────────────────
gmx editconf -f complex.gro -o complex_c.gro -c -box 3.0 3.0 3.0
gmx solvate -cp complex_c.gro -cs /usr/local/share/gromacs/top/spc216.gro \
            -box 3.0 3.0 3.0 -o solv.gro -p complex.top

# ── 4. 加离子 ───────────────────────────────────────────
gmx grompp -f $SKILL/assets/mdp/ions.mdp -c solv.gro -p complex.top -o ions.tpr
printf 'SOL\n' | gmx genion -s ions.tpr -o solv_ions.gro -p complex.top \
                 -pname NA -nname CL -neutral < /dev/null

# ── 5. 平衡（EM → NVT → NPT）─────────────────────────────
gmx grompp -f $SKILL/assets/mdp/em.mdp  -c solv_ions.gro -p complex.top -o em.tpr
gmx mdrun -s em.tpr  -deffnm em  -ntmpi 1 -ntomp 6 -pin on
gmx grompp -f $SKILL/assets/mdp/nvt.mdp -c em.gro -r solv_ions.gro -p complex.top -o nvt.tpr
gmx mdrun -s nvt.tpr -deffnm nvt -ntmpi 1 -ntomp 6 -pin on
gmx grompp -f $SKILL/assets/mdp/npt.mdp -c nvt.gro -r solv_ions.gro -t nvt.cpt -p complex.top -o npt.tpr
gmx mdrun -s npt.tpr -deffnm npt -ntmpi 1 -ntomp 6 -pin on

# ── 6. 生产 MD（100 ns）─────────────────────────────────
gmx grompp -f $SKILL/assets/mdp/md.mdp -c npt.gro -t npt.cpt -p complex.top -o md.tpr
gmx mdrun -s md.tpr -deffnm md -ntmpi 1 -ntomp 6 -pin on -x md.xtc

# ── 7. 分析 ─────────────────────────────────────────────
bash $SKILL/scripts/analyze.sh .

# ── 8. MM/PBSA 结合自由能 ────────────────────────────────
bash $SKILL/scripts/run_mmpbsa.sh . Protein lig 1 50
```

---

## 分模块使用说明

### 1. 分子对接（`scripts/dock_ligand.sh`）

```bash
bash dock_ligand.sh <受体PDB> <配体mol2/sdf> [盒子中心X Y Z] [盒子边长]
```

- **gnina 直接吃 PDB 受体**，内部自动转 PDBQT，不需要 meeko；
- 不给盒子中心时用受体几何中心做盲对接（边长自动设 30）；
- 输出多构象 SDF（用 `$$$$` 分隔），日志末尾给出
  `affinity / intramol / CNN pose score / CNN affinity` 四列；
- **挑 pose 优先看 CNN pose score**，不是 affinity。

> ⚠️ **盒子必须对准口袋**。盒子偏了软件仍会输出构象，但 affinity 全是 `0.00`
> （配体落在蛋白外），极易被误读成"结合很弱"。详见 `references/docking.md`。

### 2. 力场选型（`references/force-fields.md`）

包含：

- 本机 GROMACS 实际内置的力场清单（实测 `ls /usr/local/share/gromacs/top/*.ff`）；
- 按体系类型分支的**选型决策树**；
- 蛋白力场 × 水模型**配对表**；
- **CHARMM36 官方强制的 10 项 mdp 参数**及原因（`vdw-modifier = force-switch` 等）；
- **GROMOS 的官方警告**（按 twin-range 截断参数化，单截断下密度可能偏离）；
- AMBER → GROMACS 的两条官方转换路径（`amb2gmx.pl` / `acpype`）与**参数合并细节**；
- gmx_MMPBSA 1.7.0 中与力场相关的全部变量（`igb` / `PBRadii` / `source_force_field` / legacy 变量）。

### 3. 配体参数化与力场合并（`scripts/setup_ligand.sh`、`scripts/merge_gaff_ff.py`）

```bash
# 一条命令：SMILES/mol2 → 3D → GAFF → 合并力场 → 复合物
bash setup_ligand.sh ligand.smi topol.top protein.gro complex
```

`merge_gaff_ff.py` 也可以单独用：

```bash
python3 merge_gaff_ff.py \
  --lig-itp lig.acpype/lig_GMX.itp \
  --ff-src   /usr/local/share/gromacs/top/amber99sb.ff \
  --out-dir  ./amber99sb_gaff.ff \
  --strip-lig-itp ./lig.itp
```

> GROMACS 的拓扑**只允许一个 `[ atomtypes ]` 段**，GAFF 配体不合并进去会报
> `Atomtype c3 not found` 或 `Syntax error`。而且 GROMACS 2025 把原子类型拆在
> `ffnonbonded.itp` 和 `atomtypes.atp` **两个文件**里，脚本两个都会写。

### 4. CHARMM36 体系

先从 MacKerell 实验室下载 GROMACS 格式的 CHARMM36，装到
`/usr/local/share/gromacs/top/charmm36.ff/`，然后：

```bash
gmx pdb2gmx -f receptor.pdb -o protein.gro -i posre.itp -p topol.top \
            -ff charmm36 -water tips3p -ignh < /dev/null
```

mdp **必须**换成 `assets/mdp/charmm36_*.mdp`（已按官方 10 项参数配好，
`grompp` 实测 0 warning）。若下载包里没有 `tips3p.itp`，
可从 `/usr/local/share/gromacs/top/charmm27.ff/tips3p.itp` 复制。

### 5. MM/PBSA（`scripts/run_mmpbsa.sh`）

```bash
bash run_mmpbsa.sh <工作目录> [受体组] [配体组] [起始帧] [结束帧]
```

可选环境变量：

| 变量 | 默认 | 说明 |
|---|---|---|
| `SRC_FF` | `amber` | `source_force_field`，**只用于来源标注与提示**，不切换力场、不改拓扑 |
| `PBRADII` | `4` | 连续介质半径集；4 = mbondi3（配 `igb=8`）；7 = charmm_radii（仅 CHARMM） |
| `IGB` | `8` | 广义 Born 模型；8 = GB-Neck2（1.7.0 起默认，5 → 8 会让旧结果不可比） |
| `IDECOMP` | `3` | 能量分解方式，1/2 = 逐残基，3/4 = 成对（资源消耗大） |
| `PRINT_RES` | `within 6` | 要打印的残基；官方推荐值 |
| `GMXBIN` | `/usr/local/bin` | GROMACS 可执行文件目录 |
| `MMPBSA_HOME` | `$HOME/micromamba` | micromamba 根目录 |

> ⚠️ **读结果千万别看错段**：`mmpbsa.out` 里 `Ligand:` 段是**配体自身的绝对能量**，
> 真正的 **ΔG_bind 在 `Delta (Complex - Receptor - Ligand):` 段的 `ΔTOTAL` 行**。
> 报告不确定度用 `Block SD` / `Block SEM`，不要用 `SD(Prop.)`/`SEM(Prop.)`
> （对 1e6 量级的绝对键合能做比例传播会放大到无意义）。

### 6. 结果分析（`scripts/analyze.sh`）

```bash
bash analyze.sh <工作目录>
```

产出 RMSD、RMSF、能量分量、二级结构、Ramachandran、自由能面等图表。

两条结合能分析路线，**优先第一条**：

1. **GROMACS 2025 原生组间能量** —— mdp 里设
   `energygrps = Protein LIG Protein_LIG`，`.edr` 里直接出现
   `Coul-SR:Protein-Ligand` 等项，`gmx energy` 按编号选即可，**精确无抵消误差**；
2. **`dit xvg_energy_compute`** —— 旧式相减法，输入必须恰好 5 列且顺序为
   `Time, LJ(SR), Disper.corr., Coulomb(SR), Coul.recip.`，
   对应 `gmx energy` 选能量项 **7, 8, 9, 10**。

---

## 力场选型速查

**默认推荐（蛋白-配体小分子）**

| 部件 | 选择 | 理由 |
|---|---|---|
| 蛋白 | `amber99sb` / `amber99sb-ildn` | 与 GAFF 配体原生配套；gmx_MMPBSA 的 ParmEd 转换链路最稳 |
| 配体 | **GAFF**（acpype `-a gaff`） | GROMACS 官方推荐，与 AMBER 蛋白/核酸力场兼容 |
| 水 | `tip3p` | 内置开箱即用 |
| 静电 | `coulombtype = PME`，`rcoulomb = 1.2` | 周期体系标配 |
| 范德华 | `vdwtype = cutoff`，`rvdw = 1.2`，`cutoff-scheme = Verlet` | 配 `DispCorr = EnerPres` |

**注意**

- 做 MM/PBSA 时**蛋白力场固定用 AMBER 系**。实测 OPLS 蛋白拓扑经 ParmEd 转 AMBER
  时键合参数出错（`BOND = 609605 kcal/mol`），某些残基能量溢出，sander 打印
  `*********`，gmx_MMPBSA 解析时直接 `ValueError` 崩掉。AMBER 系是原生转换，无此问题。
- **力场与水模型必须成套选**：蛋白用 AMBER 配 `tip3p`，CHARMM 配 `tips3p`。
- **CHARMM36 必须换一整套 mdp**（官方 10 项），不能只改 `-ff`。
- **GROMOS 不建议用于定量自由能计算**（官方手册有明确警告）。

---

## 常见问题（踩坑清单）

完整清单见 `references/gromacs-2025-api.md`、`references/mmpbsa.md`、`references/docking.md`，
这里列最高频的几条：

| 现象 | 原因 | 处理 |
|---|---|---|
| `gmx pdb2gmx` 100% CPU 卡死、不报错 | 非交互环境下缺 `-ff` 或 `-water` | 两个都要给；**不要加 `-ter`**（会进交互式末端选择，同样死循环） |
| `grompp` 报 `Too many warnings (N)` | GROMACS 2025 把警告升级为错误 | 改 mdp 消除警告，别急着加 `-maxwarn` |
| `grompp` 报 `comm-mode` 相关警告 | 周期体系用了 `comm-mode = Angular` + `comm-grps = Protein` | 改成 `Linear` + `System`（模板已改好） |
| `gmx editconf -center` 行为异常 | 2025 里 `-center` 变成**向量**参数 | `-c` 才是布尔居中 |
| `grompp -o` 报错 | 2025 里 `-o` 只接受 `.tpr`，且默认输出名取拓扑文件名前 5 字符 | 显式指定 `-o`；不同体系用不同名字 |
| `trjconv` 结果异常 | `-pbc` 与 `-fit` 必须分两步做 | 组选择次数：`-pbc` 2 次，`-center -fit` 3 次 |
| `Atomtype c3 not found` | GAFF 配体原子类型没合并 | 用 `merge_gaff_ff.py` |
| `dit ndx_add` 建的组莫名其妙是空的 | **dit 不覆盖已存在的输出文件**，会另存为 `xxx_时间戳.ndx` 且仍返回 0 | 写之前先 `rm -f` 目标文件 |
| `editconf` 直接段错误 | 索引里有超出结构原子数的组 | 受体/配体索引改用**局部编号** |
| 对接 affinity 全是 `0.00` | 盒子没对准口袋 | 核对盒子中心 |
| `could not convert string to float: '*********'` | 某残基能量溢出（常见于 OPLS 拓扑转换） | 换 AMBER 系蛋白力场 |

---

## 与 GROMACS 2025 的兼容性说明

GROMACS 2025 相比 2024 有大量命令行变更，**本 skill 的所有脚本与模板都针对 2025 校正**。
最常见的几处（详见 `references/gromacs-2025-api.md`）：

- `genbox` / `tprconv` 被移除；
- `grompp` 不再接受 `-ff`；
- 警告升级为错误；
- `editconf -center` 变成向量参数；
- `gmx hbond` 被重写，`-h` 变布尔；xpm 输出要用 `gmx hbond-legacy -hbm`。

**若你要用 GROMACS 2024 或更早版本**，本 skill 的 mdp 模板需要回退
`comm-mode`/`comm-grps` 等设置，并重新验证零警告。

---

## 致谢与许可证

本项目基于以下开源项目构建，**在此致谢**：

| 项目 | 用途 | 许可证 |
|---|---|---|
| [GROMACS](https://www.gromacs.org/) | MD 主引擎 | GPL-3.0 |
| [DuIvyTools](https://github.com/CharlesHahn/DuIvyTools) | `dit` 命令；**本 skill 的 mdp 模板由其模板派生并按 GROMACS 2025 修正** | GPL-3.0 |
| [acpype](https://github.com/alanxw/acpype) | GAFF/GAFF2 配体参数化 | GPL-2.0 |
| [gmx_MMPBSA](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA) | MM/PBSA 结合自由能 | GPL-3.0 |
| [gnina](https://github.com/gnina/gnina) | CNN 打分对接 | Apache-2.0 |
| [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) | 分子对接 | Apache-2.0 |
| [Open Babel](https://openbabel.org/) | 格式转换 | MPL-2.0 |
| [meeko](https://github.com/forlilab/Meeko) | 受体 PDBQT 制备 | LGPL-3.0 |

由于 `assets/mdp/` 下的模板派生自 **GPL-3.0** 的 DuIvyTools，
本项目整体以 **GNU GPL v3** 发布。

```
Copyright (C) 2026 centuryrain1980

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```

---

## 引用与复现

如果你在论文里使用了本 skill 的流程或模板，请**同时引用**：

- GROMACS：Abraham et al., *GROMACS 2025*, J. Chem. Theory Comput. (2025)
- DuIvyTools：<https://github.com/CharlesHahn/DuIvyTools>
- acpype：Silva & Vranken, *ACPYPE*, J. Chem. Inf. Model. 55, 2467 (2015)
- gmx_MMPBSA：Valdes-Tresanco et al., *gmx_MMPBSA*, J. Comput. Chem. 36, 1483 (2015)
- gnina：McNutt et al., *gnina*, J. Chem. Inf. Model. 57, 2062 (2017)
- Vina：Trott & Olson, *AutoDock Vina*, J. Comput. Chem. 31, 455 (2010)
- 所用场：对应文献（如 AMBER ff99SB、CHARMM36 等），见 `references/force-fields.md`