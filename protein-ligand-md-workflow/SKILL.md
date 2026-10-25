---
name: protein-ligand-md-workflow
description: 蛋白-配体复合物的分子动力学全流程——分子对接、体系建模与溶剂化、力场选型（AMBER/CHARMM36/OPLS 与水模型配对）、GAFF 配体参数化、EM/NVT/NPT 平衡、生产 MD、MM/PBSA 结合自由能与结果分析（基于 GROMACS 2025 + DuIvyTools/dit + gmx_MMPBSA + gnina）。当任务涉及 vina/gnina 对接、pdb2gmx、acpype、solvate、genion、grompp、mdrun、轨迹后处理、力场与水模型选择，或需要分析 RMSD/RMSF/结合能/二级结构/Ramachandran/自由能面时使用。也用于排查 GROMACS 2025 版命令参数变更导致的失败。
---

# 蛋白-配体 MD 全流程

## 适用范围

从 PDB 结构和配体 SMILES/结构出发，完成：建模 → 溶剂化 → 参数化 → 平衡 → 生产 MD → 结果分析，产出可复现的轨迹与分析图表。

## 本机环境（已实测确认）

| 组件 | 位置 | 说明 |
|---|---|---|
| GROMACS | `/usr/local/bin/gmx` | 2025.4，AVX2_256 优化版（比 apt 版快 13.8%~19.1%） |
| GROMACS MPI | `/usr/local/bin/gmx_mpi` | 多进程版，`mpirun -np N gmx_mpi mdrun` |
| DuIvyTools | `dit` 命令 | 0.6.0，xvg/xpm 可视化与后处理 |
| 配体工具 | `~/bin/obabel`、`~/bin/vina`、`~/bin/gnina` | 格式转换、对接 |
| 配体参数化 | `acpype` | GAFF 力场，见 `scripts/setup_env.sh` |
| MM/PBSA | `gmx_MMPBSA` 1.7.0 | micromamba 环境 `mmpbsa`（Python 3.12 + AmberTools 26.0） |
| 对接 | `gnina` **1.3.3**（旧版 1.1 留档在 `~/bin/gnina_v1.1`）/ `vina` 1.2.7 | `~/bin`；gnina 直接吃 PDB 受体，推荐用 |

运行环境是 **WSL2 的 Ubuntu 26.04**。

**三条硬性环境约束**（踩过坑）：

1. **以普通用户运行**，不要用 root。root 建的目录会让后续普通用户写入失败（`Permission denied`）。
2. **`/tmp` 在 WSL 重启后会被清空**。工作目录放 `~/` 或工作区，不要放 `/tmp`。
3. **这台 CPU 是混合架构**（11 物理核 / 22 线程）。单体系推荐 `-ntomp 6 -pin on`；大体系用 `mpirun -np 4 gmx_mpi mdrun`。

## 工作流总览

每个阶段的详细命令与注意事项见对应 reference，按需读取，不要一次全读。

| 阶段 | 做什么 | 详见 |
|---|---|---|
| 0 环境与自检 | 装依赖、验证 gmx/dit 可用 | `scripts/install_all.sh`（一键装齐）、`scripts/setup_env.sh`（自检） |
| **0.5 力场选型** | **定蛋白力场 / 配体力场 / 水模型 / mdp 组合** | **`references/force-fields.md`** |
| 1 蛋白准备 | `pdb2gmx` 生成蛋白拓扑 | `references/system-setup.md` |
| 2 分子对接 | Vina / gnina 把配体对进口袋 | `references/docking.md` |
| 3 配体参数化 | GAFF/AC 配体电荷与拓扑、力场合并 | `references/system-setup.md` |
| 4 组装与溶剂化 | 复合物合并、溶剂化、加离子、索引组 | `references/system-setup.md` |
| 5 平衡 | EM → NVT → NPT | `references/equilibration.md` |
| 6 生产 MD | 长程轨迹 | `references/equilibration.md` |
| 7 结果分析 | RMSD/能量/二级结构/FEL/图表 | `references/analysis.md` |
| 7.5 结果报告 | 把分析结果汇总成带目录的 PDF（封面+指标表+配体2D+每图一页） | `scripts/make_report.py` |
| 8 结合自由能 | MM/PBSA (gmx_MMPBSA) | `references/mmpbsa.md` |

**GROMACS 2025 的命令参数有大量变更**（`editconf -center` 变成向量、`grompp` 去掉 `-ff`、警告变错误、`genbox`/`tprconv` 被移除……）。凡是在 GROMACS 2025 上报参数错误或命令卡死，**先读 `references/gromacs-2025-api.md`**，不要照搬旧教程。

## 最小可运行骨架

> **默认力场**（蛋白-配体小分子体系，见 `references/force-fields.md`）：
> 蛋白 `amber99sb` + 配体 **GAFF** + 水 `tip3p`。
> 用 AMBER 系是因为它与 GAFF 配体原生配套，且 MM/PBSA 转换链路已被 gmx_MMPBSA 官方支持。
> CHARMM36 必须改用 `assets/mdp/charmm36_*.mdp`（官方要求 10 项参数全改）。

```bash
cd ~/projects/<体系名> && mkdir -p mdp && cd mdp

# 1. 蛋白
gmx pdb2gmx -f protein.pdb -o protein.gro -i posre.itp -p topol.top \
            -ff amber99sb -water tip3p -ignh < /dev/null

# 2. 居中 + 建盒（GROMACS 2025 新语法）
gmx editconf -f protein.gro -o protein_c.gro -c -box 3.0 3.0 3.0

# 3. 溶剂化（-p 会自动更新 topol.top 的 [ molecules ]）
#    TIP3P 几何用 spc216.gro（GROMACS 不带 tip3p.gro，两者几何相同）
gmx solvate -cp protein_c.gro -cs /usr/local/share/gromacs/top/spc216.gro \
            -box 3.0 3.0 3.0 -o solv.gro -p topol.top

# 4. 加离子
gmx grompp -f mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr
printf 'SOL\n' | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
                 -pname NA -nname CL -neutral < /dev/null

# 5. 平衡 + 生产（模板在 assets/mdp/，已按 GROMACS 2025 修正）
gmx grompp -f mdp/em.mdp  -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -s em.tpr  -deffnm em  -ntmpi 1 -ntomp 6 -pin on
gmx grompp -f mdp/nvt.mdp -c em.gro  -r solv_ions.gro -p topol.top -o nvt.tpr
gmx mdrun -s nvt.tpr -deffnm nvt -ntmpi 1 -ntomp 6 -pin on
gmx grompp -f mdp/npt.mdp -c nvt.gro -r solv_ions.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -s npt.tpr -deffnm npt -ntmpi 1 -ntomp 6 -pin on
gmx grompp -f mdp/md.mdp  -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -s md.tpr  -deffnm md  -ntmpi 1 -ntomp 6 -pin on -x md.xtc
```

## 必须遵守的硬性约束

- **`grompp` 不能有警告**。GROMACS 2025 把警告升级为错误（`Too many warnings (N)` 直接终止）。真遇到无害警告时用 `-maxwarn N`，但优先改 mdp 消除它。
- **不要用 `comm-mode = Angular` + `comm-grps = Protein`**。周期体系下会触发 2 个警告（旋转移除有伪影、原子未归属任何组）导致 grompp 失败。用 `comm-mode = Linear` + `comm-grps = System`（assets 里的模板已改好）。
- **`.xtc` 需要 mdp 里显式写 `nstxout-compressed`**，只写 `nstxout` 只会得到 `.trr`。
- **`grompp -o` 只接受 `.tpr`**，且**输出文件名默认取拓扑文件名前 5 个字符**（`topol_big.top` → `topol.tpr`）。不同体系务必显式 `-o` 指定，避免互相覆盖。
- **配体参数化用 acpype 2026 语法**：`acpype -i lig.mol2 -b lig -a gaff -o gmx -c bcc`。
  旧写法（`-o lig.gro -n lig`）会报错；且 `-n` 现在是净电荷，不是分子名。
- **GAFF 配体 + OPLS 蛋白必须合并力场**：用 `scripts/merge_gaff_ff.py`，
  否则 GROMACS 会报 `Syntax error`（重复的 atomtypes 段）或 `Atomtype c3 not found`。
- **`dit ndx_add -c a-b` 的区间是左闭右开**（等同 Python `range(a,b)`）。要包含原子 b 必须写 `a-(b+1)`。
  这个坑与 DuIvyTools 自己的文档不一致，务必核对生成出的组实际原子数。
- **MM/PBSA 的三份索引必须用局部编号**：受体索引写 `1..N_prot`、配体索引写 `1..N_lig`。
  用复合物编号会让 GROMACS 2025 的 `editconf` 直接**段错误**（不是报错）。
- **`dit ndx_add` 不会覆盖已存在的输出文件**：它会另存为 `xxx_时间戳.ndx` 并且**仍返回 0**，
  新建的组就这样静默丢失（后面报 `Define a valid index group`）。
  写之前必须先 `rm -f` 目标文件。
- **组名查找要大小写不敏感**：GROMACS 本身不区分大小写（`Non-Water` == `non-Water`），
  但 shell 里用 awk 精确匹配 `[ LIG ]` 时，索引里若是 `[ lig ]` 就会数出 0 个原子。
- **MM/PBSA 的 ΔG 在输出 `Delta (Complex - Receptor - Ligand):` 段的 `ΔTOTAL` 行**，
  `Ligand:` 段是配体自身的绝对能量（不是 ΔG）。详见 `references/mmpbsa.md`。
- **对接盒子必须对准口袋**，否则 gnina/Vina 仍会输出构象但 affinity 全是 0.00（配体落在蛋白外）。
- **交互式选择必须喂输入**。`gmx energy` / `trjconv` / `rms` / `genion` 等会提示选择，用 `printf 'Protein\nProtein\n' | gmx ...` 或 `< /dev/null` 喂空行。
- **给 `grompp` 传了 `-n index.ndx` 就等于完全替换默认索引组**，自定义文件里必须包含 mdp 引用的每一个组名。
- **`gmx pdb2gmx` 在非交互环境下必须同时给 `-ff` 和 `-water`**，否则会静默死循环（100% CPU、不报错、不结束）。
  同理**不要加 `-ter`**：它会进入交互式末端选择，喂 `< /dev/null` 同样死循环。
- **gnina v1.3.3 启动需要 CUDA 运行时库**（`libcudnn.so.9` 等）。
  已用 pip 装到 `~/.local/lib/python3.14/site-packages/nvidia/`，
  由包装脚本 `~/apps/gnina/gnina` 自动设 `LD_LIBRARY_PATH`。
  直接调 `gnina_v1.3.3` 原始二进制会报 `error while loading shared libraries: libcudnn.so.9`。
- **不同 gnina 版本的 CNN 分数不能横向比较**（v1.3 重训并蒸馏了 CNN 模型，数值尺度变了）。
  同一批配体要统一版本；详见 `references/docking.md`。
- **力场与水模型必须成套选**：蛋白用 AMBER 系就配 `tip3p`，CHARMM 系配 `tips3p`，
  水模型跨力场混用会让静电参数与力场假设不自洽。`gmx pdb2gmx -ff` 的可选值就是
  `/usr/local/share/gromacs/top/*.ff` 的目录名，本机清单见 `references/force-fields.md` 第 1 节。
- **CHARMM36 必须用 `assets/mdp/charmm36_*.mdp`**：官方要求 `vdw-modifier = force-switch`、
  `rvdw-switch = 1.0`、`DispCorr = no` 等 10 项，通用模板不满足。
- **MacKerell 的 `charmm36.ff` 要先打段名补丁**：`ethers.n.tdb` 里的 `[ MET1 ]` / `[ HYD1 ]`
  会被 GROMACS 误当成"给蛋白 MET/HYD 用的末端补丁"，导致
  `atom C1 not found in buiding block 1MET while combining tdb and rtp`。
  改名成 `[ ETHERS_MET1 ]` / `[ ETHERS_HYD1 ]` 即可（`install_all.sh --with-charmm36` 会自动打）。
  注意**必须改成不以 MET/HYD 开头**的名字。详见 `references/force-fields.md` 第 4.2 节。
- **CHARMM36 力场目录只能有一处**：GROMACS 找力场的顺序是"当前目录 → `GMXLIB` → 安装树"。
  放 `~/apps/` 下必须设 `GMXLIB`（否则报 `Could not find force field 'charmm36'`）；
  但**系统力场库和 `GMXLIB` 不能同时有**同名力场，否则报
  `Force field 'charmm36' occurs in 2 places, but not in the current directory.`。
  推荐 `sudo cp` 到 `/usr/local/share/gromacs/top/`，任何场景都不用管环境变量。
- **GROMOS 不建议用于定量自由能计算**：官方手册警告它按 twin-range 截断参数化，
  在 GROMACS 的单截断下密度等物理性质可能偏离。
- **`micromamba create` 会先清空目标环境再下载**：镜像一旦 404 就留下一个空壳环境
  （`envs/mmpbsa/bin/` 变成 0 个文件），原有环境直接报废。
  所以判断"环境是否存在"要看**关键可执行文件**（`python` / `gmx_MMPBSA` / `cpptraj`），
  **不能**只看 `micromamba env list` 的名字——它的输出行首有两个空格，
  `grep -qE "^mmpbsa"` 匹配不到，会误判成"不存在"从而触发重建。
  `scripts/setup_mmpbsa.sh` 已按此修正。
- **生产 MD 被中断后用 `--resume` 续跑**：`run_simulation.sh <dir> <nsteps> --resume` 会用
  `-cpi md.cpt -append` 从检查点继续并**追加**到原轨迹（不会覆盖）。
  想彻底重跑用 `--restart`（会先删掉旧的 `md.*` 产物）。
  不带这两个参数时默认从头跑，会覆盖旧轨迹。
- **复合物拓扑 `complex.top` 里写的是合并力场目录的绝对路径**，
  所以**整个工作目录要整体保留/拷贝**，只拷 `complex.top` 到别处会导致
  `grompp` 找不到力场文件。
- **MM/PBSA 的 `mmpbsa.in` 不要写 `forcefields` / `protein_forcefield` /
  `ligand_forcefield` / `ions_parameters`**：gmx_MMPBSA 1.7.0 官方明确这些是 legacy 变量，
  新文件不要加；`source_force_field` 只用于来源标注，不会切换力场。

## 分析环节

DuIvyTools 的 `dit` 命令负责 xvg/xpm 的可视化与后处理，`scripts/analyze.sh` 是可复用的分析脚本。
蛋白-配体结合能分解有两条路，优先用第一条：

1. **GROMACS 2025 原生组间能量**——mdp 里设 `energygrps = Protein Ligand Protein_Ligand`，
   .edr 里直接出现 `Coul-SR:Protein-Ligand` 等项，`gmx energy` 按编号选即可，精确无抵消误差。
2. **`dit xvg_energy_compute`**——旧式相减法（prolig − pro − lig），
   输入必须是恰好 5 列且顺序为 Time, LJ(SR), Disper.corr., Coulomb(SR), Coul.recip.，
   对应 `gmx energy` 选能量项 **7, 8, 9, 10**。

详见 `references/analysis.md`。

## 参考文件

- `references/gromacs-2025-api.md` —— **命令参数变更与排错**，遇到报错先读
- `references/system-setup.md` —— 蛋白/配体准备、组装、溶剂化、离子、索引组
- `references/equilibration.md` —— EM/NVT/NPT/生产 MD 的参数原理与调优
- `references/analysis.md` —— 各类分析的 GROMACS 命令 + `dit` 命令映射表
- `references/force-fields.md` —— **力场选型与参数转换**：本机可用力场清单、决策树、水模型配对、
  CHARMM36 官方强制 mdp、AMBER→GROMACS 转换细节、gmx_MMPBSA 力场相关参数、常见错误速查
- `references/mmpbsa.md` —— MM/PBSA 结合自由能：环境、命令行参数硬性要求、踩过的坑
- `references/docking.md` —— 分子对接 (Vina/gnina)：盒子定位、PDBQT 制备、坑
- `assets/mdp/` —— 已修正的 em/ions/nvt/npt/md 参数模板（AMBER/OPLS 通用）
- `assets/mdp/charmm36_*.mdp` —— CHARMM36 专用模板（官方 10 项强制参数已配好）
- `scripts/install_all.sh` —— **一键安装全部运行环境**（系统包 → GROMACS 2025 编译 → Python 工具链 → MM/PBSA 环境 → 对接工具 → PATH 持久化 → 自检），幂等，新机器首选
- `scripts/mdflow.py` —— **统一命令行入口**：`verify / prepare / solvate / simulate / analyze / mmpbsa / report / all`，支持 `--dry-run` 预演
- `scripts/make_report.py` —— 生成 PDF 报告（figurebook）；无 reportlab 时自动降级为 HTML
- `scripts/setup_env.sh` —— 环境自检与依赖安装
- `scripts/merge_gaff_ff.py` —— 把 GAFF 配体原子类型并入蛋白力场（GROMACS 只允许一个 atomtypes 段，必须合并）
- `scripts/assemble_complex.py` —— 蛋白 + 配体 → 复合物坐标与拓扑（含残基号/原子序号重编）
- `scripts/dock_ligand.sh` —— 分子对接 (gnina + CNN 重打分)
- `scripts/setup_ligand.sh` —— 配体路径一条命令：SMILES → 3D → GAFF → 力场合并 → 复合物
- `scripts/setup_mmpbsa.sh` —— 安装 gmx_MMPBSA 环境（Python 3.12 + AmberTools）
- `scripts/run_mmpbsa.sh` —— 一键跑 MM/PBSA 并汇总结果
- `scripts/run_simulation.sh` —— 批量执行平衡+生产阶段
- `scripts/analyze.sh` —— 批量分析并出图

DuIvyTools 本身是 GPLv3 开源项目（https://github.com/CharlesHahn/DuIvyTools），本 skill 的 mdp 模板在其基础上按 GROMACS 2025 修正。
