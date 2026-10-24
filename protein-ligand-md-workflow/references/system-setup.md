# 体系构建：蛋白、配体、组装、溶剂化

## 目录

- [1. 蛋白准备](#1-蛋白准备)
- [2. 配体准备与参数化](#2-配体准备与参数化)
- [3. 复合物组装](#3-复合物组装)
- [4. 溶剂化](#4-溶剂化)
- [5. 加离子](#5-加离子)
- [6. 索引组](#6-索引组)
- [7. 拓扑合并要点](#7-拓扑合并要点)

---

## 1. 蛋白准备

用 `pdb2gmx` 把 PDB 转成 GROMACS 拓扑。

```bash
gmx pdb2gmx -f protein.pdb -o protein.gro -i posre.itp -p topol.top \
            -ff oplsaa -water spc -ignh < /dev/null
```

| 选项 | 作用 | 备注 |
|---|---|---|
| `-ff` | 力场 | `oplsaa` / `amber99sb-ildn` / `charmm36-jul2021` |
| `-water` | 水模型 | 必须与后续 solvate 用的模板一致（`spc` 对应 `spc216.gro`） |
| `-ignh` | 忽略输入氢 | PDB 一般无氢，保留 |
| `-i posre.itp` | 位置限制 | **NVT/NPT 阶段要用 `-r` 传入对应坐标** |
| `-p topol.top` | 拓扑 | 会自动写入 `[ molecules ]` 段 |

产出：

- `protein.gro` —— 加氢后的坐标
- `topol.top` —— 蛋白拓扑（后续要追加配体、溶剂）
- `posre.itp` —— 位置限制参数

**常见问题**

- 缺主链原子 → `atom C used in that entry is not found in the input file`。检查 PDB 是否被裁剪过。
- 非标准残基 → 报 `residue XXX not found in residue topology database`。
  要么改成标准残基，要么用 `gmx pdb2gmx -res` 配合自定义残基库。
- 多构象/altLoc → 先用 `pdb2pdb` 或编辑 PDB 选定单一构象。

**已有蛋白拓扑时**（比如从 PDBBind 取的 AMBER 力场体系）可跳过本步，直接进第 3 步。

---

## 2. 配体准备与参数化

这是蛋白-配体流程里最容易出错的一环。本节命令均在**本机实测通过**（以阿司匹林为测试配体）。

### 2.1 从 SMILES 得到 3D 结构

**Open Babel 3.x 的语法变了**：`-o` 只接输出**格式**，文件名要用大写 `-O`。

```bash
# 错误(2.x 写法): obabel: cannot write output format!
obabel -i smi -o mol.sdf --gen3d lig.smi

# 正确(3.x 写法):
obabel -ismi lig.smi -osdf -Omol.sdf --gen3d
# 或先转 mol2(acpype 只吃 pdb/mdl/mol2/SMILES, 不吃 sdf)
obabel -ismi lig.smi -omol2 -Olig.mol2 --gen3d
```

本机 `obabel` 来自 pip 的 `openbabel-wheel`，程序在 `/usr/local/bin/obabel`；
`~/bin/obabel` 是一个转发脚本。**若 `~/bin/obabel` 报找不到 `apps/obabel/usr/bin/obabel`，
说明该路径已失效，改用 `/usr/local/bin/obabel`。**

### 2.2 质子化与手性

```bash
obabel -imol2 -omol2 -Olig_pH.sdf --pH 7.4
```

有手性中心时**必须人工确认**质子化态与立体化学。常见问题：咪唑环（His 类似）、
胍基（Arg）、羧酸（Asp/Glu）在 pH 7.4 的电离态。

### 2.3 GAFF 参数（acpype 2026）

**acpype 2026 重构了命令行**，旧写法会报错：

```bash
# 错误(2023 写法): Invalid value for '-o': 'lig.gro' is not one of 'all','gmx','cns','charmm'
acpype -i mol.sdf -o lig.gro -a gaff -n lig

# 正确(2026 写法):
acpype -i lig.mol2 -b lig -a gaff -o gmx -c bcc
```

| 选项 | 2026 语义 | 备注 |
|---|---|---|
| `-i` | 输入文件（.pdb/.mdl/.mol2）或 SMILES 字符串 | **不收 .sdf** |
| `-b` | 输出文件名前缀（basename） | 取代旧的 `-o 文件名` |
| `-a` | 原子类型：`gaff` / `gaff2` / `amber` / `amber2` | 默认是 **gaff2** |
| `-c` | 电荷方法：`gas`/`bcc`/`abcg2`/`user` | 默认 bcc |
| `-n` | **净电荷（整数）** | 旧版这里填的是分子名！ |
| `-o` | 输出拓扑格式：`all`/`gmx`/`cns`/`charmm` | 不再是文件名 |

产物在 `<basename>.acpype/` 目录：

| 文件 | 内容 |
|---|---|
| `lig_GMX.itp` | GROMACS 拓扑，**自带 `[ atomtypes ]` 段**（GAFF 类型） |
| `lig_GMX.gro` | GROMACS 坐标 |
| `lig_GMX_OPLS.itp` | 变体，类型名换成 `opls_XXX`（见 2.5 的坑） |
| `lig_bcc_gaff.mol2` | 带 GAFF 类型与 BCC 电荷的 mol2 |
| `posre_lig.itp` | 配体位置限制 |

安装（Ubuntu 26.04 的 `rich` 包由 Debian 管理，直接 pip 装会冲突）：

```bash
pip3 install --break-system-packages --ignore-installed rich \
    -i https://pypi.tuna.tsinghua.edu.cn/simple acpype mdanalysis mdtraj
```

装好后 `acpype` 在 `~/.local/bin`（不在默认 PATH，需 `export PATH="$HOME/.local/bin:$PATH"`）。

### 2.4 检查参数质量

```bash
grep -A9 "\[ atoms \]" lig.acpype/lig_GMX.itp | head -9
```

关注：原子类型是否被正确识别（不应出现 `C1`/`C2` 这类通用碳）、电荷总和是否等于净电荷。
实测阿司匹林：9 个 GAFF 类型（c3 c o os ca oh hc ha ho），BCC 电荷总和 **-0.0000**。

### 2.5 GAFF 原子类型并入蛋白力场（关键难点）

**问题**：GROMACS 的拓扑只允许一个 `[ atomtypes ]` 段。acpype 的 `lig_GMX.itp` 自带一段，
直接 `#include` 会报：

```
Fatal error:
Syntax error - File lig.itp, line 3
```

**GROMACS 2025 的力场目录把原子类型写在两个文件里**，必须同时改，只改一个会报
`Atomtype c3 not found`：

| 文件 | 格式 |
|---|---|
| `ffnonbonded.itp` | `[ atomtypes ]` 段：`name bond_type mass charge ptype sigma epsilon` |
| `atomtypes.atp` | `name mass ;` 列表 |

**解法**：用 skill 提供的脚本做合并（已实测）：

```bash
python3 <skill>/scripts/merge_gaff_ff.py \
    --lig-itp  lig.acpype/lig_GMX.itp \
    --ff-src   /usr/local/share/gromacs/top/oplsaa.ff \
    --out-dir  ./oplsaa_gaff.ff \
    --strip-lig-itp ./lig.itp
```

它会：复制力场目录 → 把 GAFF 类型追加到 `ffnonbonded.itp` 和 `atomtypes.atp` →
生成去掉了 `[ atomtypes ]` 段的 `lig.itp`。

然后拓扑里指向**绝对路径**的合并力场：

```
#include "/绝对路径/oplsaa_gaff.ff/forcefield.itp"
#include "lig.itp"
```

> **不要用 acpype 的 `lig_GMX_OPLS.itp` 变体**：它把 GAFF 类型映射成 `opls_XXX`，
> 但遇到 OPLS 没有对应的 GAFF 类型时会生成占位符 `opls_x`（质量 0），报
> `Atomtype opls_x not found`。实测阿司匹林的酯基氧就踩了这个坑。

**验证**：合并后跑一次 `gmx grompp` + 最小化，若报 `Atomtype ... not found` 说明
`atomtypes.atp` 没同步改；若报 `[ bondtypes ]` 语法错误，说明 `ffnonbonded.itp`
被改坏（用 `scripts/merge_gaff_ff.py` 重做即可）。

## 3. 复合物组装

把蛋白和配体放进同一个坐标文件，并合并拓扑。

### 3.1 坐标合并

蛋白和配体要在**同一坐标系**。若配体来自对接结果（Vina/gnina），坐标已对齐，直接合并。

```bash
# 用 gmx concatenate 或手工 cat（注意 gro 头部原子数与末行盒矢量）
gmx editconf -f protein.gro -o prot_only.gro
```

手工合并 gro 时要改三处：标题行、原子数行、末行盒矢量。

更稳妥的方式是让配体**最后加入**，然后重新生成盒：

```bash
gmx editconf -f complex.gro -o complex_c.gro -c -box 3.0 3.0 3.0
```

### 3.2 拓扑合并

在 `topol.top` 里插入配体的 `.itp`，并更新 `[ molecules ]`：

```
#include "lig.itp"

...

[ molecules ]
; Compound        #mols
Protein_chain_A     1
LIG                 1
SOL               442
```

配体在 `[ molecules ]` 里的**顺序必须和坐标文件里的原子顺序一致**。
蛋白在前、配体在后，是最不容易出错的安排。

---

## 4. 溶剂化

```bash
gmx solvate -cp complex_c.gro -cs /usr/local/share/gromacs/top/spc216.gro \
            -box 3.0 3.0 3.0 -o solv.gro -p topol.top
```

| 选项 | 说明 |
|---|---|
| `-cp` | 溶质（蛋白+配体复合物） |
| `-cs` | 溶剂模板，水模型要与 pdb2gmx 的 `-water` 一致 |
| `-box X Y Z` | 盒子尺寸（nm） |
| `-p topol.top` | **自动更新** `[ molecules ]` 里的溶剂数量 |
| `-scale` | 溶剂模板缩放，默认 0.57 |

**盒子尺寸经验值**：蛋白/配体复合物距离盒边缘至少留 1.0–1.2 nm（大于 vdW 截断半径），
否则会有严重的周期性相互作用伪影。溶质最大线性尺寸 + 2.4 nm 是安全下界。

**水的密度**：solvate 填出来的水密度是标准值（约 33.4 个水分子/nm³）。

---

## 5. 加离子

```bash
# 用 ions.mdp 建 tpr
gmx grompp -f mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr

# 按电荷中和加离子
printf 'SOL\n' | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
                 -pname NA -nname CL -neutral < /dev/null
```

`genion` 会**替换**水分子而不是新增，所以总原子数不变。

| 选项 | 说明 |
|---|---|
| `-neutral` | 中和体系净电荷 |
| `-conc` | 指定浓度（mol/L），如 `-conc 0.15` |
| `-pname/-nname` | 正/负离子名 |
| `-rmin` | 离子与溶质的最小距离，默认 0.6 nm |

离子命名必须和拓扑里的 `[ moleculetype ]` 对应（oplsaa 用 `NA`/`CL`；CHARMM 用 `NA`/`CL`；
AMBER 体系常见 `Na`/`Cl` 或 `K`/`Cl`）。genion 会自动改写 `topol.top` 的 `[ molecules ]`。

**加离子后必须重新最小化**——否则带电粒子落在不合适的位置会有巨大斥力。

---

## 6. 索引组

### 生成默认组

```bash
printf 'q\n' | gmx make_ndx -f solv_ions.gro -o index.ndx
```

默认组已包含 `Protein`、`Water`、`non-Water`、`Water_and_ions` 等，**通常够用**，不必手工维护。

### 需要自定义组时

用 `dit`（DuIvyTools）比手写文件可靠：

```bash
# 新增组：-al 组名，-c 原子序号（支持 1-10-3 这样的步进写法）
dit ndx_add -f index.ndx -al LIG -c 1301-1340 -o index_new.ndx

# 把一个组均分成 N 份
dit ndx_split -f index.ndx -al Protein 2 -o index_split.ndx

# 查看有哪些组
dit ndx_show -f index.ndx
```

### 常用组名约定

mdp 里引用的组名必须存在。推荐统一命名：

| 组名 | 含义 | 用途 |
|---|---|---|
| `Protein` | 蛋白（默认组） | RMSD、对齐 |
| `LIG` | 配体 | 配体 RMSD、接触分析 |
| `Water` | 水 | tc-grps |
| `non-Water` | 非水 | tc-grps |
| `Protein_Ligand` | 蛋白+配体 | 能量分解 |

> 组名匹配**不区分大小写**，mdp 里写 `Non-Water` 也能匹配 `non-Water`。

---

## 7. 拓扑合并要点

最终的 `topol.top` 结构：

```
#include "/usr/local/share/gromacs/top/oplsaa.ff/forcefield.itp"
#include "/usr/local/share/gromacs/top/oplsaa.ff/spc.itp"
#include "lig.itp"

[ system ]
Protein-ligand complex in water

[ molecules ]
; Compound        #mols
Protein_chain_A     1
LIG                 1
SOL               442
NA                  4
CL                  2
```

**一致性检查清单**

- `[ molecules ]` 里的分子种类与顺序 == 坐标文件里的分子顺序
- 每个分子的数量 == 坐标文件里该分子的实际数量
- 配体的 `[ moleculetype ]` 名（acpype 的 `-n`）== `[ molecules ]` 里写的名字

校验方法：

```bash
# 坐标文件的分子数（标题行）
head -1 solv_ions.gro
# grompp 会做原子数匹配，数目不符时报：
#   number of coordinates in coordinate file (x.gro, N) does not match topology (topol.top, M)
```
