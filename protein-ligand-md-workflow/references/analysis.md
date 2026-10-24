# 结果分析

分析分两层：**GROMACS 自带工具**负责算，**DuIvyTools（`dit`）**负责画图和后处理。

---

## 目录

- [1. 分析前的准备](#1-分析前的准备)
- [2. 常用分析项](#2-常用分析项)
- [3. 结合能分解](#3-结合能分解)
- [4. 二级结构与 Ramachandran](#4-二级结构与-ramachandran)
- [5. 自由能面 (FEL)](#5-自由能面-fel)
- [6. dit 命令速查表](#6-dit-命令速查表)
- [7. 绘图与出图](#7-绘图与出图)
- [8. 分析脚本](#8-分析脚本)

---

## 1. 分析前的准备

### 轨迹预处理

去掉周期性边界伪影——蛋白分子会跳到盒子另一侧，直接算 RMSD 会出现巨大尖峰。

GROMACS 2025 要求 **PBC 处理与旋转拟合分两步执行**，合并会报
`PBC condition treatment does not work together with rotational fit`。

```bash
# 第一步: 处理周期性边界（需要 2 次组选择）
printf 'Protein\nProtein\n' | gmx trjconv -s md.tpr -f md.xtc \
    -o whole.xtc -pbc whole

# 第二步: 居中 + 旋转拟合（需要 3 次组选择）
printf 'Protein\nProtein\nProtein\n' | gmx trjconv -s md.tpr -f whole.xtc \
    -o traj_fit.xtc -center -fit rot+trans
```

### 选择分析组

所有分析命令都要选原子组，脚本化时统一喂：

```bash
printf 'Protein\nProtein\n' | gmx rms -s md.tpr -f md.xtc -o rms.xvg
```

---

## 2. 常用分析项

| 分析 | 命令 | 产出 |
|---|---|---|
| 均方根位移 | `gmx rms -s md.tpr -f md.xtc -o rms.xvg -tu ps` | 稳定性 |
| 残基 RMSF | `gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res` | 柔性 |
| 回转半径 | `gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg` | 紧密程度 |
| 溶剂可及面积 | `gmx sasa -s md.tpr -f md.xtc -o sasa.xvg` | 埋藏面积 |
| 氢键 | `gmx hbond -s md.tpr -f md.xtc -num hbond_num.xvg` | 氢键统计 |
| 氢键图谱 | `gmx hbond-legacy -s md.tpr -f md.xtc -hbm hbond.xpm` | 氢键存在性矩阵（2025 版 xpm 要用 legacy） |
| 接触数 | `gmx mindist` / `gmx distance` | 配体-残基距离 |
| 二面角 | `gmx angle`, `gmx dihedral` | 特定几何量 |
| 聚类 | `gmx cluster -s md.tpr -f md.xtc -cl clusters.xpm -g cluster.log -dist rmsd-dist.xvg` | 构象多样性 |
| 主成分分析 | `gmx covar -s md.tpr -f md.xtc -o eigenval.xvg -v eigenvec.trr -av avg.pdb` | 构象空间 |
| 熵/焓 | `gmx gangle`/`gmx enthalpy` | 热力学量 |
| 自由能面 | `gmx gnmats` + `gmx xpm` | 见第 5 节 |

**配体 RMSD**（判断配体是否稳住在口袋里）：

```bash
# 需要先把配体单独抽出来作为拟合参考
printf 'LIG\nLIG\n' | gmx rms -s md.tpr -f md.xtc -o rms_lig.xvg -fit rot+trans
```

---

## 3. 结合能分解

有两条路，**优先用第一条**。

### 路线 A：GROMACS 2025 原生组间能量（推荐）

mdp 里设置能量组：

```
energygrps             = Protein LIG Protein_LIG
```

跑完 MD 后，.edr 里直接含有组间能量项：

```
48 Coul-SR:Protein-Protein    49 LJ-SR:Protein-Protein
52 Coul-SR:Protein-rest       53 LJ-SR:Protein-rest
56 Coul-SR:rest-rest          57 LJ-SR:rest-rest
```

在蛋白-配体体系里，`Protein_LIG` 与 `rest` 的组合会给出
`Coul-SR:Protein_LIG-rest` 这类项，直接选出来即可：

```bash
# 先看有哪些组间能量项
printf 'q\n' | gmx energy -f md.edr -o /dev/null 2>&1 | grep -E ':'

# 选编号导出（编号见 gromacs-2025-api.md 的完整列表）
printf '52\n53\n' | gmx energy -f md.edr -o inter.xvg
```

优点：**没有大数相减的抵消误差**，数值直接来自力计算。

### 路线 B：`dit xvg_energy_compute`（相减法）

```bash
# 三个文件都必须恰好 5 列，顺序为 Time, LJ(SR), Disper.corr., Coulomb(SR), Coul.recip.
# 对应 gmx energy 的能量项 7, 8, 9, 10
printf '7\n8\n9\n10\n' | gmx energy -f md.edr -o prolig.xvg -g g1.gro
printf '7\n8\n9\n10\n' | gmx energy -f md.edr -o pro.xvg    -g g2.gro
printf '7\n8\n9\n10\n' | gmx energy -f md.edr -o lig.xvg    -g g3.gro

dit xvg_energy_compute -f prolig.xvg pro.xvg lig.xvg -o binding.xvg
```

注意 `gmx energy` 的 `-g` 选项在 2025 版已被移除，上面的 `-g` 会报
`Unknown command-line option -g`，**去掉 `-g` 即可**。

输出 `binding.xvg` 的列：

```
Time, LJ(SR), Disper.corr., Coulomb(SR), Coul.recip.,
LJ(all), Coulomb(all), Short-Range(all), Long-Range(all), Total Energy
```

其中 `Total Energy` 就是结合能的负值（prolig − pro − lig）。

> `dit xvg_energy_compute` 会严格校验列名与列数，legend 顺序不对会报
> `the legend order should be: Time, LJ(SR), Disper.corr., Coulomb(SR), Coul.recip.`

### MM/PBSA 结合自由能

本机已装好完整环境（micromamba 环境 `mmpbsa`：Python 3.12 + AmberTools 26.0 + gmx_MMPBSA 1.7.0）。

```bash
export MAMBA_ROOT_PREFIX=$HOME/micromamba
$HOME/micromamba/bin/micromamba run -n mmpbsa gmx_MMPBSA \
  -cs complex.pdb -ci idx_c.ndx -cg "Protein" "lig" -ct com.xtc -cp complex.top \
  -rs rec.pdb -ri idx_r.ndx -rg "Protein" -rt rec.xtc -rp rec.top \
  -ls lig.pdb -li idx_l.ndx -lg "lig" -lt lig.xtc -lp lig.top \
  -i mmpbsa.in -o mmpbsa.out
```

一键脚本：`bash scripts/run_mmpbsa.sh <工作目录>`。

**参数硬性要求**（每条都实测过）：结构给 `.pdb`/`.tpr`、拓扑给 `.top`、
三个部分必须用**三份不同的索引文件**且受体/配体索引用**局部编号**、
`-cg` 要给两个组名。细节见 `references/mmpbsa.md`。

---

## 4. 二级结构与 Ramachandran

### 二级结构

GROMACS 2023+ 的 `gmx dssp` 输出的是**符号矩阵**（每帧一行残基二级结构符号），
不是老的 xpm。需要 `dit dssp` 转成标准格式：

```bash
gmx dssp -s md.tpr -f md.xtc -o dssp.dat
dit dssp -f dssp.dat -o dssp_out
```

产出：

- `dssp_out.xpm` —— 二级结构随时间的二维图
- `dssp_out_sc.xvg` —— 各二级结构元件的残基编号随时间变化
- `dssp_out_residue.xvg` —— 各帧各残基的二级结构

画图：

```bash
dit xpm_show -f dssp_out.xpm -o dssp.png
```

### Ramachandran 图

```bash
printf 'Protein\nProtein\n' | gmx rama -s md.tpr -f md.xtc -o rama.xvg
dit xvg_rama -f rama.xvg -o rama.png
```

输出会给出落在一般/甘氨酸/脯氨酸区域的残基统计。

---

## 5. 自由能面 (FEL)

FEL 用两个反应坐标（常用 RMSD 与 RMSF，或两个 RMSD）做二维直方图。

```bash
# 1) 算两个反应坐标
printf 'Protein\nProtein\n' | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ps
printf 'Protein\nProtein\n' | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res -fit

# 2) gnmats 生成自由能网格
gmx gnmats -f rmsd.xvg -f2 rmsf.xvg -o gnmats.xpm -hist rmsd.xvg rmsf.xvg \
            -xmin 0 -xmax 5 -ymin 0 -ymax 5 -bins 50
# 注意 gnmats 在 2025 版可能以子命令形式提供，先用 gmx gnmats -h 确认

# 3) 画图
dit xpm_show -f gnmats.xpm -o fel.png
```

`dit xpm_show` 支持 2D 热图等高线，参数 `-xlabel/-ylabel/-zlabel`、`-xmin/-xmax` 等可调。

---

## 6. dit 命令速查表

### XVG 类（数据文件）

| 命令 | 用途 |
|---|---|
| `dit xvg_show -f x.xvg -o out.png` | 画曲线（需 `MPLBACKEND=Agg` 无头环境） |
| `dit xvg_compare -f a.xvg b.xvg` | 多条曲线对比 |
| `dit xvg_ave -f x.xvg -o ave.xvg` | 求平均 |
| `dit xvg_ave_bar -f x.xvg` | 平均值柱状图 + 误差棒 |
| `dit xvg_box_compare -f x.xvg` | 箱线图/小提琴图对比 |
| `dit xvg_show_distribution -f x.xvg` | 分布直方图 |
| `dit xvg_show_scatter -f x.xvg` | 散点图 |
| `dit xvg_show_stack -f x.xvg` | 堆叠面积图 |
| `dit xvg_combine -f a.xvg b.xvg -o c.xvg` | 合并列 |
| `dit xvg_energy_compute -f prolig.xvg pro.xvg lig.xvg` | 结合能分解 |
| `dit xvg_rama -f rama.xvg` | Ramachandran 图 |

### XPM 类（二维矩阵）

| 命令 | 用途 |
|---|---|
| `dit xpm_show -f x.xpm -o out.png` | 画二维图 |
| `dit xpm2csv -f x.xpm -o out.csv` | 转 CSV (x,y,z) |
| `dit xpm2dat -f x.xpm -o out.dat` | 转 DAT (N*N) |
| `dit xpm_diff -f a.xpm b.xpm` | 两图差值 |
| `dit xpm_merge -f a.xpm b.xpm` | 两图各取一半合并 |

### 其他

| 命令 | 用途 |
|---|---|
| `dit mdp_gen` | 生成 mdp 模板 |
| `dit show_style` | 查看绘图样式文件 |
| `dit find_center -f x.gro index.ndx` | 求原子组几何中心 |
| `dit dccm_ascii -f dccm.dat -o dccm.xpm` | DCCM 转 xpm |
| `dit dssp -f dssp.dat -o out` | GROMACS 2023+ 的 dssp 符号矩阵转 xpm/xvg |
| `dit ndx_add -f in.ndx -al NAME -c 1-10,20-30 -o out.ndx` | 新增索引组 |
| `dit ndx_split -f in.ndx -al NAME 2 -o out.ndx` | 拆分索引组 |
| `dit ndx_show -f in.ndx` | 列出所有组 |

常用参数：`-eg {matplotlib,plotext,plotly,gnuplot}` 选绘图引擎；
`-x/-y/-z` 坐标轴标签；`-t` 标题；`-b/-e` 起止；`-csv` 导出 CSV。

---

## 7. 绘图与出图

### 无头环境

```bash
export MPLBACKEND=Agg     # 无显示器时必须，否则 matplotlib 会报错或挂起
```

### 字体与中文

`dit` 的绘图样式文件在 `DuIvyTools/data/mplstyle/`。中文标签可能缺字体，
建议图内用英文，注释用中文。

### 分辨率与尺寸

```bash
dit xvg_show -f rms.xvg -o rms.png --figsize 8 6
```

---

## 8. 分析脚本

`scripts/analyze.sh` 封装了常用分析（轨迹预处理 → RMSD/RMSF/回转半径/SASA/氢键 →
二级结构 → Ramachandran → 能量分解 → 出图），可直接运行：

```bash
bash scripts/analyze.sh <工作目录>
```

其中 `<工作目录>` 需包含 `md.tpr`、`md.xtc`、`md.edr`（或 `topol.tpr`/`md.xtc` 等）。
