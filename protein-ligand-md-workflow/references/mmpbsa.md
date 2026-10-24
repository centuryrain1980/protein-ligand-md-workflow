# MM/PBSA 结合自由能 (gmx_MMPBSA)

## 目录

- [环境](#环境)
- [命令行参数硬性要求](#命令行参数硬性要求)
- [输入文件 mmpbsa.in](#输入文件-mmpbsain)
- [完整流程](#完整流程)
- [踩过的坑](#踩过的坑)
- [结果解读](#结果解读)

---

## 环境

> **选型提示**：做 MM/PBSA 时**蛋白力场固定用 AMBER 系**（`amber99sb` / `amber99sb-ildn`），
> 配体用 GAFF。实测 OPLS 蛋白拓扑经 ParmEd 转 AMBER 时键合参数出错
> （BOND = 609605 kcal/mol），某些残基能量溢出导致 sander 输出 `*********`、gmx_MMPBSA 解析崩溃。
> 完整依据见 `references/force-fields.md`。

gmx_MMPBSA 是独立项目（不在 GROMACS 源码树里，2025 版已移除 `src/contrib/gmx_MMPBSA`）：

```bash
git clone https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git
```

**硬性依赖**（缺一不可）：

| 依赖 | 要求 | 说明 |
|---|---|---|
| Python | **>=3.11, <3.13** | 本机系统是 3.14，装不了；必须用独立环境 |
| AmberTools | 提供 `cpptraj` / `tleap` / `sander` / `parmchk2` | conda-forge 的 `ambertools` 包（286 MB） |
| parmed | >=4.2.2 | pip 可装 |
| GROMACS | 在 PATH 中 | 通过 `gmx_path` 指定 |

**本机已装好的环境**：

```bash
export MAMBA_ROOT_PREFIX=$HOME/micromamba
$HOME/micromamba/bin/micromamba run -n mmpbsa gmx_MMPBSA -h
```

- Python 3.12.15 / cpptraj V7.6.2 (AmberTools 26.0) / tleap / sander / parmed 4.3.1
- gmx_MMPBSA 1.7.0
- 重建环境：`bash scripts/setup_mmpbsa.sh`

---

## 命令行参数硬性要求

这是最容易卡住的地方。以下每一条都实测确认：

| 参数 | 格式要求 |
|---|---|
| `-cs` / `-rs` / `-ls`（结构） | 只接受 **`.pdb` 或 `.tpr`**，给 `.gro` 报 `does not correspond to the required structure format ['.tpr', '.pdb']` |
| `-cp` / `-rp` / `-lp`（拓扑） | 只接受 **`.top`**，给 `.tpr` 报 `does not correspond to the required structure format ['.top']` |
| `-ci` / `-ri` / `-li`（索引） | **必须是三个不同的文件**，内容相同也会报 `Several args are duplicated in the command-line` |
| `-cg`（复合物组） | **需要两个组名**：受体组 + 配体组，如 `-cg Protein lig`。只给一个会吞掉后面的参数（argparse 的 list 行为） |
| `-rg` / `-lg` | 各一个组名 |
| `-ct` / `-rt` / `-lt`（轨迹） | **必须与对应结构的原子数一致**。复合物轨迹只能含复合物原子，受体轨迹只能含受体原子 |

### 索引必须是"局部编号"

这是第二个大坑：`-ri` / `-li` 用的索引里的原子号必须是**相对于该结构自身**的编号，
不能直接用复合物编号。工具内部会执行：

```
gmx editconf -f lig.pdb -n <配体索引> -o _GMXMMPBSA_LIG.pdb
```

若索引里写的是复合物编号（如 1232-1252），而 `lig.pdb` 只有 21 个原子，
**GROMACS 2025 的 editconf 会直接段错误（Segmentation fault），而不是报错**。

正确做法：受体索引写 `1..N_prot`，配体索引写 `1..N_lig`。

---

## 输入文件 mmpbsa.in

Fortran Namelist 格式，分 `&general` / `&gb` / `&decomp` 三段。

```fortran
&general
  sys_name      = "SystemName",
  startframe    = 1,
  endframe      = 50,
  interval      = 1,
  source_force_field = "amber",
  PBRadii       = 4,
  temperature   = 298.15,
  gmx_path      = "/usr/local/bin",
  solvated_trajectory = 0,
  verbose       = 2,
/

&gb
  igb         = 8,
  saltcon     = 0.150,
  intdiel     = 1.0,
  extdiel     = 78.5,
  surften     = 0.0072,
/

&decomp
  idecomp     = 3,
  print_res   = "within 6",
/
```

要点：

- **`idecomp` 不能为 0**（当 `-cg` 给了两个组时即启用了分解），报
  `IDECOMP cannot be 0 for Decomposition analysis!`
- `print_res` 格式是 `链/残基`，如 `"R/1-75"`、`"A/2-10,35,41 B/104"`；
  不带链前缀会报 `We expected something like this: A/2-10,35,41 B/104`
- 分解必须同时给出受体和配体残基，否则报
  `For decomposition analysis, you must define residues for both receptor and ligand!`
- `igb=8` 是 **GB-Neck2**（`igb=1` 才是 HCT）；`igb=2/5` 是 GB-OBC1/OBC2，更快但精度低一些

### 1.7.0 的默认值变了（写方法学时必须记下来）

| 变量 | 旧默认 | **1.7.0 默认** | 影响 |
|---|---|---|---|
| `igb` | 5 | **8** | GB-OBC2 → GB-Neck2。旧脚本不显式写 `igb`，结果会变 |
| `PBRadii` | 3 | **4** | mbondi2 → mbondi3，与 `igb=8` 配套 |
| `pb exdi` | 80.0 | **78.5** | 溶剂介电常数 |

**报告结果时必须写明实际用的 `igb` 和 `PBRadii`**，否则无法复现。

`igb` 与 `PBRadii` 的官方配对表：

| `igb` | 模型 | 配套 `PBRadii` |
|---|---|---|
| 1 | GB-HCT | 2 (mbondi) |
| 2 | GB-OBC1 | 3 (mbondi2) |
| 5 | GB-OBC2 | 3 (mbondi2) |
| 7 | GB-Neck | 1 (bondi) |
| **8** | **GB-Neck2** | **4 (mbondi3)** ← 默认 |

CHARMM 体系另用 `PBRadii = 7`（`charmm_radii`，官方注明只配 CHARMM 力场制备的体系）。

### 不要写的 legacy 变量

官方原文：`forcefields`、`protein_forcefield`、`ligand_forcefield`、`ions_parameters`
是 legacy 的拓扑制备设置，**只为兼容旧输入文件保留，不参与当前的拓扑式计算，
"Do not add them to new input files"**。

同时，旧的 `gmx_MMPBSA -lm ligand.mol2`（用 tleap 现算配体参数）**这条路已经移除**：
GROMACS 计算路径必须用 `-cp` 提供复合物拓扑，小分子配体必须已经在拓扑里。

### `source_force_field` 的真实语义

官方原文：*"used for provenance and advisories … **never selects a radius set or changes a topology**"*。
它只是告诉程序"这是 AMBER/CHARMM/OPLS 系的体系"，好给出对应提醒，
**不能**靠改它来切换力场（力场在 MD 拓扑里就定死了）。取值：
`auto`(默认) / `amber` / `charmm` / `opls` / `gromos` / `other`。

### `idecomp` 与 `print_res` 的选择

| `idecomp` | 模式 |
|---|---|
| 1 | 逐残基分解，1-4 项计入内部能 |
| 2 | 逐残基分解，1-4 EEL/VDW 分别计入 EEL/VDW |
| 3 | **成对**分解，1-4 项计入内部能 |
| 4 | **成对**分解，1-4 EEL/VDW 分别计入 EEL/VDW |

`print_res` 默认 `"within 6"`（受体与配体间距 6 Å 内的残基），官方推荐用它。
官方明确 `print_res = "all"` *"generally not recommended"*（多数残基贡献为 0，浪费时间）。

官方警告：`idecomp = 3/4`（成对模式）在大体系 + 大范围 `print_res` 时会产生
**巨大的临时 mdout**（官方举例：250 个残基约 500 MB，因为是 62500 对），
解析时内存占用也很大。

力场/半径集/水模型的完整选型依据见 `references/force-fields.md`。

生成模板：

```bash
gmx_MMPBSA --create_input gb      # 生成 mmpbsa.in
gmx_MMPBSA --input-file-help       # 列出所有可用变量
```

---

## 完整流程

```bash
# 1. 拆受体/配体结构与轨迹
gmx editconf -f complex.gro -o complex.pdb
gmx editconf -f rec.gro    -o rec.pdb
gmx editconf -f lig.gro    -o lig.pdb
printf "Protein\n" | gmx trjconv -s md.tpr -f md.xtc -n idx_c.ndx -o rec.xtc
printf "lig\n"     | gmx trjconv -s md.tpr -f md.xtc -n idx_c.ndx -o lig.xtc
printf "Protein_lig\n" | gmx trjconv -s md.tpr -f md.xtc -n idx_c.ndx -o com.xtc

# 2. 三份索引(复合物用复合物编号, 受体/配体用局部编号)

# 3. 受体/配体各自的 .top(只列自己的分子)

# 4. 运行
gmx_MMPBSA \
  -cs complex.pdb -ci idx_c.ndx -cg "Protein" "lig" -ct com.xtc -cp complex.top \
  -rs rec.pdb    -ri idx_r.ndx -rg "Protein"    -rt rec.xtc -rp rec.top \
  -ls lig.pdb    -li idx_l.ndx -lg "lig"        -lt lig.xtc -lp lig.top \
  -i mmpbsa.in -o mmpbsa.out
```

一键脚本：`bash scripts/run_mmpbsa.sh <工作目录> [受体组] [配体组] [起始帧] [结束帧]`

---

## 踩过的坑

| 现象 | 原因 | 处理 |
|---|---|---|
| `Define a valid index group`（MMPBSA_Error） | 索引里的受体/配体组是**空的** | 见下面两条：①`dit ndx_add` **不会覆盖已存在的输出文件**——它会另存为 `xxx_时间戳.ndx` 并且**仍返回 0**，新组就这样静默丢了；写之前必须先 `rm -f` 目标文件。②组名查找用了**大小写敏感**的 awk，而 GROMACS 本身**不区分大小写**（索引里是 `[ lig ]`、参数传 `LIG` 就查不到，会数出 0 个原子） |
| 配体组原子数 = 0 | 同上 | 同上；`run_mmpbsa.sh` 已改为每次重建索引 + 大小写不敏感查找，并在组为空时直接报错 |
| `Several args are duplicated in the command-line` | 三个部分用了同一份索引/轨迹/拓扑文件 | 拆成三份独立文件 |
| `does not correspond to the required structure format ['.tpr','.pdb']` | `-cs` 给了 `.gro` | 转成 `.pdb` |
| `does not correspond to the required structure format ['.top']` | `-cp` 给了 `.tpr` | 给 `.top` |
| `error: unrecognized arguments: md.xtc` | `-cg` 只给一个组，把后面的文件名吃掉了 | `-cg` 给两个组名 |
| `Group Protein referenced in the .mdp file was not found` | mdp 里的组名不在索引中 | 补组 |
| **Segmentation fault**（editconf） | 索引里有超出结构原子数的组 | 受体/配体索引改用局部编号 |
| `IDECOMP cannot be 0` | 给了两个组即启用分解 | `idecomp = 1..4` |
| `We expected something like this: A/2-10,35,41 B/104` | `print_res` 缺链前缀 | `"R/1-75"` |
| `could not convert string to float: '*********'` | sander 分解表某个残基能量溢出 Fortran 字段宽度 | 见下 |

### `*********` 溢出

sander 的分解表首列是残基的总能量，若某残基能量超过 1e8 会打印 `*********`，
gmx_MMPBSA 解析时直接抛 `ValueError`。

常见原因：**OPLS 蛋白拓扑被 ParmEd 转成 AMBER 时键合参数出错**，导致某些残基内部能量爆炸
（实测 OPLS 体系 BOND = 609605 kcal/mol；LYS76 直接溢出）。

**解决：用 AMBER 系力场建蛋白**（`amber99sb`），转换是原生的，键合能正常：

```bash
gmx pdb2gmx -f protein.pdb -o protein.gro -p topol.top \
            -ff amber99sb -water tip3p -ignh < /dev/null
```

注意 TIP3P 的几何用 `spc216.gro` 提供（GROMACS 不带 `tip3p.gro`，两者几何相同）。

另外 `idecomp = 1`（逐残基分解）比 `idecomp = 3/4` 更容易触发该溢出，
实测 `idecomp = 3` / `4` 可以正常跑完。

---

## gmx_MMPBSA 1.7.0 的输出长什么样（实测，最容易读错的地方）

`mmpbsa.out` 在 `GENERALIZED BORN:` 之后有 **4 个小节**，每节一张表：

```
Complex:    Energy Component  Average  SD(Prop.)  SD  SEM(Prop.)  SEM  Block SD  Block SEM
            BOND / ANGLE / DIHED / VDWAALS / EEL / 1-4 VDW / 1-4 EEL / EGB / ESURF
            GGAS / GSOLV / TOTAL
Receptor:   同上
Ligand:     同上
Delta (Complex - Receptor - Ligand):   同上，但每一行前面带 Δ 前缀(ΔBOND/ΔTOTAL ...)
```

> ### ⚠️ **`Ligand:` 段是"配体自身的绝对能量"，不是 ΔG！**
>
> 真正的 **ΔG_bind 在 `Delta (Complex - Receptor - Ligand):` 段的 `ΔTOTAL` 行**。
>
> 验算（本机实测数据）：
>
> ```
> ΔBOND   = 1385219.45 − 1385229.80 − 4.26       = −14.61   (表里 ΔBOND = −14.60) ✓
> ΔEEL    =   −5471.83 −  −5512.54 − 43.89      =   −3.18   (表里 ΔEEL  =  −3.19) ✓
> ΔEGB    =   −1228.95 −  −1220.93 − (−17.86)   =    9.84   (表里 ΔEGB  =   9.85) ✓
> ΔESURF  =      35.83 −     33.92 −   2.22     =   −0.31   (表里 ΔESURF=  −0.31) ✓
> ΔTOTAL  =  1382957.53 − 1383028.72 − (−61.81)  =   −9.38   (表里 ΔTOTAL=  −9.38) ✓
> ```
>
> 也就是说：**`Ligand:` 段的 TOTAL（这里是 −61.81）是"配体自己有多稳定"，
> 而 ΔG_bind（这里是 −9.38）才是"结合有多强"。**
> 把前者当成 ΔG 汇报会得到大一个量级的错误结果。

---

## 结果解读

- **ΔG_binding = TOTAL(Complex) − TOTAL(Receptor) − TOTAL(Ligand)**，
  在输出里就是 `Delta (Complex - Receptor - Ligand):` 段的 `ΔTOTAL` 行。
- 绝对能量（BOND 等）可以非常大（本例 1.4e6 kcal/mol），但差值是有意义的。
- `GGAS` = 除 EGB 之外的所有项；`GSOLV` = EGB + ESURF；`TOTAL` = GGAS + GSOLV。
- **不要看 `SD(Prop.)` / `SEM(Prop.)`**：它们是对 1e6 量级的绝对键合能做比例传播，
  会放大到毫无意义（本例 ΔBOND 的 SEM(Prop.) 高达 2.9e5）。
  **报告不确定度用 `Block SD` / `Block SEM`**（官方也推荐 block 估计）。
- **每一帧的 ΔG 都有误差（kT 量级），需要多帧平均**，通常取 50–200 帧。
  帧数太少时 Block SEM 会比均值还大（本例 6 帧就是如此），**这样的结果不能作为结论**。

实测（泛素 + 阿司匹林，**6 帧**，igb=8，amber99sb）：

```
ΔG_bind  =  -9.38 kcal/mol   (Block SD 30.27, Block SEM 21.40)
ΔBOND    = -14.60    ΔVDWAALS = -0.91    ΔEEL  = -3.19
ΔEGB     =   9.85    ΔESURF   = -0.31
```

> 这个数只能当作"脚本跑通了"的验证，**不能作为结合自由能的结论**：
> 6 帧的 Block SEM (21.40) 比均值 (−9.38) 还大，说明采样远不够。
> 正式计算请用 50–200 帧重新跑。
