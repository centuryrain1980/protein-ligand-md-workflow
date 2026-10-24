# 力场选型与参数转换指南

> **本文件的依据**（不是经验之谈，都是官方文档 + 本机实测）：
>
> 1. GROMACS 官方手册《Force fields in GROMACS》
>    https://manual.gromacs.org/current/user-guide/force-fields.html
> 2. gmx_MMPBSA 1.7.0 官方文档（输入文件变量与兼容性说明）
>    https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/
> 3. AmberTools 官方文档（`antechamber` / `amb2gmx.pl` / `acpype`）
> 4. 本机实测：`ls /usr/local/share/gromacs/top/*.ff`、`gmx pdb2gmx -h`、
>    `~/mmpbsa_src/gmx_MMPBSA-master/docs/input_file.md`
>
> 术语白话：**力场（force field）**= 一套描述原子间作用力的参数表（键长、键角、二面角、
> 电荷、范德华半径等）。**力场文件（.ff 目录）**= GROMACS 里按力场名组织的一堆 `.itp` 参数文件。

---

## 目录

- [0. 结论先行：一套可以直接抄的配置](#0-结论先行一套可以直接抄的配置)
- [1. 本机 GROMACS 2025.4 实际内置的力场清单](#1-本机-gromacs-20254-实际内置的力场清单)
- [2. 选型决策树](#2-选型决策树)
- [3. 蛋白力场 × 水模型 配对表](#3-蛋白力场--水模型-配对表)
- [4. 各力场必须配套的 mdp 参数](#4-各力场必须配套的-mdp-参数)
  - [4.3 CHARMM27 与 CHARMM36 到底差在哪](#43-charmm27-与-charmm36-到底差在哪本机逐项实测)
- [5. 参数转换细节（AMBER → GROMACS）](#5-参数转换细节amber--gromacs)
- [6. gmx_MMPBSA 1.7.0 里与力场有关的参数](#6-gmx_mmpbsa-170-里与力场有关的参数)
- [7. 推荐配方与完整命令](#7-推荐配方与完整命令)
- [8. 常见错误速查表](#8-常见错误速查表)

---

## 0. 结论先行：一套可以直接抄的配置

**蛋白 + 小分子配体（本 skill 的主场景）默认用这一套：**

| 部件 | 选择 | 为什么 |
|---|---|---|
| 蛋白 | `amber99sb`（保守、兼容性最好）或 `amber99sb-ildn`（更新） | AMBER 系与 GAFF 配体是**原生配套**；GROMACS → ParmEd(AMBER) 的转换链路已被 gmx_MMPBSA 官方支持，本机实测无残基能量溢出 |
| 配体（小分子） | **GAFF**（acpype 默认）或 GAFF2 | GROMACS 官方明确："GAFF … 参数适合小分子，且与 AMBER 蛋白/核酸力场兼容" |
| 水 | `tip3p` | 内置、开箱即用、与 AMBER 配套最久 |
| 静电 | `coulombtype = PME`，`rcoulomb = 1.2` | 周期边界体系标配 |
| 范德华 | `vdwtype = cutoff`，`rvdw = 1.2`，`cutoff-scheme = Verlet` | 配 `DispCorr = EnerPres` |
| 键合约束 | `constraints = h-bonds`，`constraint_algorithm = lincs` | 允许 2 fs 步长 |
| 温度耦合 | `tc-grps = Water Non-Water`，`tau_t = 0.1 0.1` | 分组耦合比全体耦合更准 |

**⚠️ 特别提醒（实测踩过）**：OPLS 蛋白拓扑经 ParmEd 转成 AMBER 时**键合参数会出错**
（实测 BOND = 609605 kcal/mol），某些残基的能量超出 Fortran 字段宽度，
sander 输出 `*********`，gmx_MMPBSA 解析时直接 `ValueError` 崩掉。
**AMBER 系是原生转换，无此问题。** 所以本 skill 的 MM/PBSA 路径固定用 AMBER 系。

**⚠️ 用 CHARMM36 必须换一整套 mdp**：官方规定 10 项参数（`vdw-modifier = force-switch`、
`rvdw-switch = 1.0`、`DispCorr = no` 等），少改一项结果就不对。
本 skill 已提供 `assets/mdp/charmm36_*.mdp`（见第 4 节）。

**❌ 定量自由能计算不要用 GROMOS**：官方手册有明确警告，见第 4.3 节。

---

## 1. 本机 GROMACS 2025.4 实际内置的力场清单

`gmx pdb2gmx -ff <名字>` 的 `<名字>` 就是 `/usr/local/share/gromacs/top/` 下的
`<名字>.ff` 目录名。本机实测结果：

| `-ff` 取值 | 对应文献版本 | 该目录内置的水模型 | 评价 |
|---|---|---|---|
| `amber99sb` | AMBER ff99SB | spc, spce, tip3p, tip4p, tip4pew, tip5p | ✅ **蛋白首选** |
| `amber99sb-ildn` | ff14SB / ff19SB 的 ILDN 修正版 | 同上 | ✅ 新体系首选 |
| `amber99` | AMBER 99 | 同上 | ⚠️ 太老，隐式配体参数不全 |
| `amber96` / `amber94` / `amber03` | 老 AMBER | 同上 | ❌ 一般不用 |
| `amberGS` | AMBER GS | 同上 | ⚠️ 特殊用途 |
| `charmm27` | CHARMM27（官方移植并支持） | spc, spce, tip3p, tip4p, **tips3p**, tip5p | ✅ 官方支持；水用 `tips3p`（CHARMM modified TIP3P） |
| `oplsaa` | OPLS-AA / AA-M | spc, spce, tip3p, tip4p, tip4pew, tip5p, tip5pe | ⚠️ 小分子友好，但 MM/PBSA 转换有坑 |
| `gromos43a1` `43a2` `45a3` `53a5` `53a6` `54a7` | GROMOS 系列 | spc, spce, tip3p, tip4p | ⚠️ 官方警告，见 4.3 |

### 本机**没有**、需要额外安装的

| 缺什么 | 为什么要 | 怎么装 |
|---|---|---|
| `charmm36` | 目前主流的 CHARMM 蛋白力场 | 从 MacKerell 实验室网站下载 GROMACS 格式版本，放到 `/usr/local/share/gromacs/top/charmm36.ff/`（见 4.2 节） |
| `opc` / `opc3` | AMBER19SB 官方指定搭配的水模型 | GROMACS 未内置，需第三方文件；本机用 `tip3p` 代替 |
| `amber19sb` 独立目录 | — | GROMACS 2025.4 **没有**独立的 `amber19sb.ff`；带 ILDN 修正的 ff14SB/ff19SB 由 `amber99sb-ildn` 覆盖 |

> **关于 AMBER14SB / AMBER19SB**：GROMACS 官方手册把 AMBER94/96/99/99SB/99SB-ILDN/03/14SB/19SB/GS
> 都列为"原生支持"，但 2025.4 实际安装目录里只有 `amber03 / amber94 / amber96 / amber99 /
> amber99sb / amber99sb-ildn / amberGS` 七个。**用 `-ff amber19sb` 会直接报找不到目录**。
>
> 另外官方同一句写道：*"If you intend to use AMBER19SB in GROMACS, please use the OPC or OPC3
> model for water molecules."*（用 AMBER19SB 就该配 OPC/OPC3 水）——本机没有 OPC，
> 所以实际落地用 `amber99sb-ildn` + `tip3p`，属于常见的工程折中，**论文里要写清楚用的什么**。

---

## 2. 选型决策树

```
你的体系是什么？
│
├─ 蛋白 + 小分子配体（GAFF/GAFF2）        ← 本 skill 主场景
│   └─► amber99sb-ildn 或 amber99sb  +  GAFF 配体  +  tip3p
│       └─ 需要 MM/PBSA？ → 强制 amber99sb + tip3p（见第 6 节）
│
├─ 蛋白 + 核酸
│   └─ amber99sb-ildn（DNA/RNA 参数齐全）
│
├─ 蛋白 + 脂膜
│   └─ amber99sb-ildn（lipid 需 CHARMM36 lipid 目录，见 4.2）
│
├─ 小分子药物（无蛋白）
│   └─ oplsaa（OPLS-AA 对小分子参数最全）或 charmm27
│
├─ CHARMM36 是团队标准 / 要和 CHARMM-GUI 体系对照
│   └─ 装 charmm36.ff + 用 assets/mdp/charmm36_*.mdp（10 项参数必须全改）
│
└─ 复现某篇论文
    └─ 严格照抄论文的力场 + 水 + mdp 参数，不要"换个更好的"
```

**一条通用原则**：**同一个体系里蛋白、核酸、配体、水必须是一套自洽的搭配**，
不要蛋白用 AMBER、配体用 OPLS 却指望它自己兼容——GROMACS 只允许一个 `atomtypes` 段，
混用必然报错（用 `scripts/merge_gaff_ff.py` 合并，见 5.3 节）。

---

## 3. 蛋白力场 × 水模型 配对表

| 蛋白力场 | 官方/惯例配套水 | 本机 `-water` 可选值 | 说明 |
|---|---|---|---|
| `amber99sb` / `amber99sb-ildn` | **TIP3P**（AMBER 传统） | `tip3p` ✅ `spc` `spce` `tip4p` `tip4pew` `tip5p` | 首选 `tip3p`；TIP3P 几何用 `spc216.gro` 模板（GROMACS 不带 `tip3p.gro`，两者几何相同） |
| `charmm27` | **CHARMM modified TIP3P** | `tips3p` ✅（该目录专有）`tip3p` `spc` `spce` | 用 `tips3p`，它在水分子 H 上带 LJ 项，和标准 TIP3P 不是一回事 |
| `charmm36`（外部安装） | **CHARMM modified TIP3P** | 随下载包提供 | 见 4.2 |
| `oplsaa` | TIP3P / TIP4P / SPC 均可 | `tip3p` `tip4p` `spc` `spce` `tip4pew` `tip5p` `tip5pe` | OPLS 对 SPC/E 兼容性较好 |
| `gromos*` | SPC/E | `spc` `spce` | 见 4.3 的警告 |

---

## 4. 各力场必须配套的 mdp 参数

### 4.1 AMBER / OPLS（通用模板，`assets/mdp/{em,ions,nvt,npt,md}.mdp`）

本 skill 的 `assets/mdp/` 模板就是按这一套写的，已通过 GROMACS 2025「零警告」验证：

```
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.2
rvdw                    = 1.2
vdwtype                 = cutoff
DispCorr                = EnerPres
constraints             = h-bonds
constraint_algorithm    = lincs
comm-mode               = Linear
comm-grps               = System
```

> `DispCorr = EnerPres` 是 AMBER/OPLS 的标准做法（把长程色散修正加到压力上）。
> **CHARMM36 必须改成 `no`**（官方要求）。

### 4.2 CHARMM36（官方强制要求，10 项一个都不能少）

官方手册原文给出的 mdp 设置：

```
constraints = h-bonds
cutoff-scheme = Verlet
vdwtype = cutoff
vdw-modifier = force-switch
rlist = 1.2
rvdw = 1.2
rvdw-switch = 1.0
coulombtype = PME
rcoulomb = 1.2
DispCorr = no
```

要点解释：

- `vdw-modifier = force-switch` —— CHARMM 的范德华函数在 1.0–1.2 nm 之间做**平滑切换**
  （不是突变截断，也不是加校正项），这是 CHARMM 参数化时假定的前提。
  用普通的 `vdwtype = cutoff` 会得到不同的能量。
- `rvdw-switch = 1.0` + `rvdw = 1.2` —— 切换起点与截断点。
- `DispCorr = no` —— **不要**加色散修正项，因为 force-switch 本身已经处理了。
  官方补充："色散修正在脂质**单层**时需要，双层（bilayer）不需要。"
- `rlist = 1.2` —— 邻居表半径。

#### ⚠️ 必读：charmm36.ff 需要打一个段名补丁才能用（本机实测根因）

**先说结论**：CHARMM36 **可以**在 GROMACS 2025 上正常使用，
但 MacKerell 原始的 `charmm36.ff` 直接用 `pdb2gmx` 会报：

```
Fatal error:
atom C1 not found in buiding block 1MET while combining tdb and rtp
```

这不是"新 GROMACS 用不了新力场"，而是 **CHARMM36 力场文件里的一个 tdb 段名
与 GROMACS 的命名约定撞车**，打一行补丁即可解决。

##### 根因（已定位到具体文件与代码）

1. `charmm36.ff/ethers.n.tdb`（CHARMM36 特有的**醚类/糖脂**末端补丁）里有两个段：

   ```
   [ MET1 ]     ; Append terminal methyl group adjacent to CH2
   [ HYD1 ]     ; Complete terminal methyl group adjacent to ether O
   ```

   在 **CHARMM** 的语义里，`MET1`/`HYD1` 表示"比标准 MET/HYD **多一个原子**的残基变体"
   （用于醚连接的糖脂头基）。

2. 在 **GROMACS** 的 tdb 语义里，段名会被拆成"残基名 + 数字"，
   `MET1` 被理解成 **"给 MET 残基用的末端补丁"**。

3. 于是 GROMACS 把这个**醚脂质**的补丁误套到了**蛋白**的 MET 残基上
   （泛素恰好以 MET 开头，是第 1 号残基 → 报错里的 "1MET"）。
   该补丁要连接到一个叫 `C1` 的原子，而蛋白 MET 的 rtp 里根本没有 C1
   → GROMACS 在 `pdb2top.cpp` 里 `gmx_fatal` 退出。

   ```cpp
   /* src/gromacs/gmxpreprocess/pdb2top.cpp */
   if (patch->type() == MoleculePatchType::Add) {
       gmx_fatal(FARGS, "atom %s not found in buiding block %d%s "
                         "while combining tdb and rtp", ...);
   }
   ```

4. **这不是 GROMACS 2024/2025 引入的问题**：实测 GROMACS **v2023.5** 的
   `pdb2top.cpp` 这段代码与 v2025.4 **完全一致**。
   CHARMM36 在 GROMACS 生态里一直是用 **CHARMM-GUI** 建模产出的拓扑
   （绕开了 `pdb2gmx`），所以这个坑很少有人撞上。

##### 修复（一行 sed，改名而不是删文件）

```bash
sed -i 's/\[ MET1 \]/[ ETHERS_MET1 ]/; s/\[ HYD1 \]/[ ETHERS_HYD1 ]/' \
      <charmm36.ff>/ethers.n.tdb
```

改名而不是删掉整个 `ethers.n.tdb`，是为了保留 ethers 的其它段。

> ⚠️ **注意**：段名必须改成**不以 MET / HYD 开头**的名字。
> 本机实测：改成 `[ MET1_DISABLED ]` 仍然失败（因为 GROMACS 按前缀匹配），
> 改成 `[ ETHERS_MET1 ]` 才成功。

##### 修复后的实测结果（本机）

```
gmx pdb2gmx -ff charmm36 -water tip3p     exit=0  WARNING=0
gmx solvate                                exit=0
gmx grompp (charmm36_ions.mdp)             exit=0  WARNING=0
gmx genion                                 exit=0  2728 原子
gmx grompp (charmm36_em.mdp)                exit=0  WARNING=0  ERROR=0
gmx grompp (charmm36_nvt.mdp)               exit=0  WARNING=0  ERROR=0
gmx grompp (charmm36_npt.mdp)               exit=0  WARNING=0  ERROR=0
gmx grompp (charmm36_md.mdp)                exit=0  WARNING=0  ERROR=0
```

N 端电荷检查（CHARMM36 的 NH3+ 末端，正确）：

```
1  NH3  1  MET  N    -0.30   14.007
2  HC   1  MET  H1    0.33    1.008
3  HC   1  MET  H2    0.33    1.008
4  HC   1  MET  H3    0.33    1.008
5  CT1  1  MET  CA    0.21   12.011
6  HB1  1  MET  HA    0.10    1.008
```

##### 代价

改名后 **CHARMM36 的醚连接糖脂（glycolipid）末端补丁失效**。
纯蛋白 / 蛋白-配体体系完全不受影响；
确实要建糖脂体系，请走 **CHARMM-GUI**。

##### 三条路怎么选

| 方案 | 说明 |
|---|---|
| **A. `charmm36.ff` + 打补丁（推荐，已实测）** | 打完补丁后 GROMACS 2025 上完全正常，CHARMM36 蛋白 + CGenFF 小分子自洽 |
| **B. GROMACS 自带的 `charmm27.ff`** | 免打补丁，开箱可用；但蛋白参数是 CHARMM22，且**没有**配套的通用小分子力场 |
| **C. CHARMM-GUI 建模** | 社区标准做法，产出 GROMACS 拓扑；需要糖脂/特殊修饰时用这个 |

#### 如果要用 CHARMM36 的力场文件

```bash
# 1) 从 MacKerell 实验室下载(含 GROMACS 格式)
#    https://mackerell.umaryland.edu/charmm_ff.shtml
#    -> toppar_c36_*.tgz 里的 non_charmm/gromacs/charmm36-*.ff.tgz
mkdir -p ~/apps/charmm36 && cd ~/apps/charmm36
curl -fLO "https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_feb26.tgz"
tar xzf toppar_c36_feb26.tgz
tar xzf toppar/non_charmm/gromacs/charmm36-jul2022.ff.tgz
mv charmm36-jul2022.ff charmm36.ff

# 2) 【必须】打段名补丁, 否则 pdb2gmx 会报 "atom C1 not found in buiding block 1MET"
sed -i 's/\[ MET1 \]/[ ETHERS_MET1 ]/; s/\[ HYD1 \]/[ ETHERS_HYD1 ]/' charmm36.ff/ethers.n.tdb

# 3) 装到 GROMACS 力场库(需要 sudo), 或免 sudo 用 GMXLIB
sudo cp -r charmm36.ff /usr/local/share/gromacs/top/
# export GMXLIB=/usr/local/share/gromacs/top:$HOME/apps/charmm36
```

> `scripts/install_all.sh --with-charmm36` **会自动完成第 1、2 步**（含补丁）。
> `charmm36.ff` 自带 `tip3p.itp`（CHARMM modified TIP3P），不需要从别处补。
> **水模型名是 `tip3p`**；`charmm27.ff` 里才叫 `tips3p`，别混。
#### 本 skill 提供的模板（已实测）

`assets/mdp/charmm36_{em,ions,nvt,npt,md}.mdp` 已按官方 10 项参数写好，并且：

> **实测（两套力场都验过）**：GROMACS 2025.4 上
> - `charmm27.ff` + `-water tips3p`（GROMACS 自带）
> - `charmm36.ff` + `-water tip3p`（MacKerell 版，**打完段名补丁后**）
>
> 分别跑 `pdb2gmx` → `solvate` → `genion` → 对 5 个模板逐一 `gmx grompp`，
> **两套都是 exit=0、WARNING=0、ERROR=0**。

验证过程中撞到两条 GROMACS 2025 的硬性检查，**写 mdp 时必须遵守**：

| 检查 | 报错原文 | 怎么办 |
|---|---|---|
| `vdw-modifier = force-switch` 时开关点必须小于截断点 | `With switched vdw forces or potentials, rvdw-switch must be < rvdw` | `rvdw-switch = 1.0` 必须配 `rvdw = 1.2`（不能把 rvdw 写成 1.0） |
| Verlet 邻居表 + PME 下两个截断半径不能随便设 | `With Verlet lists rcoulomb!=rvdw is not supported (except for rcoulomb>rvdw with PME electrostatics)` | 设 `rcoulomb = rvdw = 1.2`（官方 CHARMM36 块本身就是 1.2/1.2） |

> 这两条在 **em / ions 阶段最容易踩**：常规模板里 EM 用的是 `rcoulomb = 1.0 / rvdw = 1.0`
> （相等所以合法），但 CHARMM36 的 force-switch 强制要求 `rvdw > rvdw-switch = 1.0`，
> 于是 `rvdw` 只能抬到 1.2，`rcoulomb` 就必须跟着抬到 1.2。
> 本 skill 的 `charmm36_em.mdp` / `charmm36_ions.mdp` 已经改好。

### 4.3 CHARMM27 与 CHARMM36 到底差在哪（本机逐项实测）

> 数据来源：本机 `/usr/local/share/gromacs/top/charmm27.ff` 与
> `~/apps/charmm36/charmm36.ff`（MacKerell `charmm36-jul2022.ff`）逐文件解析对比。
> 下方每个数字都能复现。

#### 先纠正一个常见误解

`charmm27.ff/forcefield.doc` 自己写着：

> **CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)**
> Parameters derived from **c32b1** version of CHARMM

也就是说 **"27" 指的是脂质力场的版本号，蛋白部分其实还是 CHARMM22（1998）+ CMAP（2004）**。
`charmm27.ff` ≠ "第 27 版的蛋白力场"，而是"蛋白用 CHARMM22、脂质用 CHARMM27"的打包。
`charmm36.ff/forcefield.doc` 则写明 *Parameters taken from CHARMM36 and CGenFF 4.6*。

#### 逐项对比

| 维度 | charmm27.ff | charmm36.ff |
|---|---|---|
| 蛋白参数来源 | CHARMM22 + CMAP（c32b1） | **CHARMM36**（`par_all36m_prot.prm`, jul22） |
| 原子类型总数 | 134 | **579**（450 个为其独有） |
| 蛋白 rtp 残基块数 | 45 | **541**（含各种末端/质子化变体） |
| **CMAP 骨架修正项** | **6** | **24**（4 倍） |
| 蛋白原子 **LJ（σ/ε）** | 基准 | **完全相同**（N/CA/C/O/S/H/HA… 逐个核对全一致） |
| 蛋白原子 **电荷** | 基准 | **307 个原子里 293 个完全相同**，只有**色氨酸**改了 14 个 |
| 正二面角 / 非正二面角 | 587 / 47 | 16054 / 368；**共有 458 条里 46 条（10%）参数不同** |
| **小分子** | **无** | **CGenFF 4.6**（通用小分子力场） |
| 核酸 | DNA/RNA（CHARMM22 时代参数） | 重做（`na.rtp`） |
| 脂质 | `lipids.rtp`（CHARMM27） | `lipid.rtp`（CHARMM36） |
| 额外模块 | — | 碳水、醚、金属、硅酸盐、溶剂 |
| NBFIX 存放方式 | `ffnabonded.itp` 38.6 KB + `ffnanonbonded.itp` 32.8 KB | `nbfix.itp` 16.5 KB（~111 条 `nonbond_params`） |
| 水模型名 | `tips3p` | `tip3p`（两者都是 CHARMM modified TIP3P，H 上带 LJ） |
| **GROMACS 2025 `pdb2gmx`** | ✅ 可用（实测 0 警告） | ⚠️ 需先打 `ethers.n.tdb` 段名补丁，见 4.2 节 |

#### 具体改了哪些

**电荷**——唯一有实质改动的是**色氨酸吲哚环的电荷重分配**（19 个标准残基里只改了它）：

```
TRP:  CD1  0.035→-0.150    CE2  0.130→ 0.240    CZ2 -0.115→-0.270
      CD2 -0.020→ 0.110    CE3 -0.115→-0.250    CZ3 -0.115→-0.200
      CH2 -0.115→-0.140    NE1 -0.610→-0.510
      HD1  0.115→ 0.220    HE1  0.380→ 0.370    HE3  0.115→ 0.170
      HH2  0.115→ 0.140    HZ2  0.115→ 0.160    HZ3  0.115→ 0.140
```

ALA / ARG / LYS / ASP / GLU / TYR / PHE / SER / THR / LEU / ILE / VAL / MET / CYS / ASN / GLN / GLY / PRO
——**电荷逐原子完全相同**。

**二面角**——共有 458 条正二面角里 46 条（10%）参数不同，例如：

```
CC-CT2-NH1-C :  k = 0.16736  →  8.368      (5 倍)
```

**CMAP**——6 项 → 24 项，意味着骨架 φ/ψ 的二维势能修正覆盖面大得多，
对长程二级结构比例（螺旋/无规卷曲）的影响比对结合自由能更明显。

#### 对蛋白-配体体系意味着什么

**蛋白部分的差异其实很小**：LJ 参数完全一样，电荷只有色氨酸变了，二面角 10% 不同。
**对结合自由能数值的影响大致在 1–2 kcal/mol 量级**，一般不会推翻结论。

**真正的差别在"配套"**：

- CHARMM36 自带 **CGenFF** —— 小分子有官方通用参数，与蛋白是**自洽的一套**；
- `charmm27.ff` **没有**通用小分子力场，配体只能外挂 GAFF 或 CGenFF，
  蛋白-配体的参数一致性就弱了；
- CMAP 24 vs 6 —— 做二级结构、螺旋稳定性分析时差别会被放大。

> ⚠️ **别混用**：蛋白用 `charmm27.ff` 而配体参数来自 CHARMM36 的 CGenFF，
> 或蛋白用 CHARMM36 而配体用 GAFF —— 参数不是一套自洽的东西，
> 强行组合出来的结果没有意义。

#### 选型建议

| 你的诉求 | 建议 |
|---|---|
| 要在 GROMACS 2025 上直接跑通、且不改任何东西 | `charmm27.ff` + `-water tips3p`（本 skill 已实测 0 警告） |
| 要严格的 CHARMM36 + CGenFF 配套 | `charmm36.ff` + `-water tip3p`，**先打 4.2 节的段名补丁**（本 skill 已实测全链路 0 警告） |
| 要建糖脂 / 醚连接脂质 | 用 **CHARMM-GUI** 建模产出 GROMACS 拓扑（段名补丁会让 ethers 末端补丁失效） |
| 做 MM/PBSA 结合自由能 | 优先 **AMBER99SB + GAFF**（本 skill 默认）；CHARMM 体系要用 `PBRadii = 7`（charmm_radii） |
| 复现某篇论文 | 严格照抄论文用的力场，不要"挑更新的" |
### 4.4 GROMOS（官方警告，定量计算请避开）

GROMACS 官方手册原文：

> **Warning** — The GROMOS force fields have been parametrized with a physically incorrect
> multiple-time-stepping scheme for a twin-range cut-off. When used with a single-range cut-off
> (or a correct Trotter multiple-time-stepping scheme), physical properties, such as the density,
> might differ from the intended values.

白话：GROMOS 的参数是在**双截断（twin-range）**体系下调好的，GROMACS 用的是**单截断**，
两者不等价，模拟出来的密度等物理性质可能偏离预期。

参考：[GitLab Issue 2884](https://gitlab.com/gromacs/gromacs/-/issues/2884)、
[ChemRxiv 10.26434/chemrxiv.11474583.v1](https://doi.org/10.26434/chemrxiv.11474583.v1)。

**实务建议**：GROMOS 是联合原子（没有非极性氢）力场，跑得快，但本 skill 的目标
是蛋白-配体定量结合自由能，**建议不用**。若必须用，先在同等条件下验证密度是否正常。

---

## 5. 参数转换细节（AMBER → GROMACS）

小分子配体没有 GROMACS 原生参数，必须从 AMBER 系转换过来。

### 5.1 两条官方转换路径

官方手册原文：

> There are scripts available for converting AMBER systems (set up, for example, with GAFF) to
> GROMACS (amb2gmx.pl, or ACPYPE), but they do require AmberTools installation to work.

| 工具 | 形式 | 特点 |
|---|---|---|
| `amb2gmx.pl` | AmberTools 自带 Perl 脚本 | 直接吃 tleap 建好的 AMBER 拓扑与坐标，一次转换蛋白+配体+水；输出 GROMACS 的 `.top/.itp/.gro` |
| `acpype` | Python 独立包 | 只处理**小分子配体**，输出配体自己的 `.top/.itp/.gro`，再手动并进蛋白拓扑 |

本 skill 主用 **acpype**（只需 AmberTools 的 `antechamber`，无需整套 tleap 流程）。

### 5.2 acpype 的准确用法（acpype 2026 版语法）

```bash
acpype -i ligand.mol2 \        # 输入：pdb / mdl / mol2 / SMILES
       -b ligand \             # 输出文件名前缀（不是输出文件名！）
       -a gaff \               # 力场：gaff（默认）/ gaff2 / amber / charmm ...
       -o gmx \                # 输出 GROMACS 格式
       -c bcc \                # 电荷方法：bcc（AM1-BCC）/ gas / user
       -n 0                    # ⚠️ 净电荷（不是分子名！旧版才是分子名）
```

> **acpype 版本差异是本流程最容易翻车的地方**：
> - `-i` 只接受 `pdb / mdl / mol2 / SMILES`，给别的格式会静默不产出；
> - `-b` 是**文件名**（不含扩展名），输出为 `ligand.top/.itp/.gro`；
> - `-n` 在新版本是**净电荷**，旧版本是分子名 —— 照抄旧教程会把电荷设错。
> - 默认力场是 **gaff2**，想用 GAFF 必须显式 `-a gaff`。
> - 电荷方法默认是 `gas`，推荐改 `-c bcc`（AM1-BCC 电荷质量明显更好）。

### 5.3 力场合并：GAFF 配体 + AMBER 蛋白

GROMACS 一个拓扑里**只允许一个 `[ atomtypes ]` 段**。GAFF 配体的原子类型（c1/c2/c3/ar…）
在 AMBER 蛋白力场里没有定义，直接拼会报：

```
ERROR: Atomtype c3 not found
```
或
```
Fatal error: Syntax error in file topol.top, line 5
```

解决：把 GAFF 配体的 atomtypes **合并进蛋白力场的 atomtypes 段**。
本 skill 用 `scripts/merge_gaff_ff.py` 自动完成（复制力场目录 + 写入 atomtypes + 去重）。

```bash
python3 scripts/merge_gaff_ff.py \
  --lig-itp lig.acpype/lig_GMX.itp \        # acpype 生成的配体 itp
  --ff-src   /usr/local/share/gromacs/top/amber99sb.ff \   # 蛋白力场目录
  --out-dir  ./amber99sb_gaff.ff \           # 输出：合并后的力场目录
  --strip-lig-itp ./lig.itp                 # 同时输出去掉 [ atomtypes ] 的配体 itp
```

两个必须知道的细节：

1. **GROMACS 2025 的力场目录把原子类型拆在两个文件里**，两个都要写，
   否则报 `Atomtype XXX not found`：
   - `ffnonbonded.itp` 里的 `[ atomtypes ]` 段（含 charge / ptype / sigma / epsilon）
   - `atomtypes.atp`（name + mass 列表）
2. **合并方向：以蛋白力场为主干**，把配体 atomtypes 追加进去。
   同名类型但参数不同时会被静默取"先出现的那个"，所以务必核对合并后
   `ffnonbonded.itp` 里是否有语义冲突。
3. 合并后的目录名必须以 `.ff` 结尾。`gmx pdb2gmx` 会**先在当前工作目录**找
   `<名字>.ff`，找不到再去 GROMACS 库目录找——这正是本地合并目录能生效的原理。
   `assemble_complex.py` 会把拓扑里的 `#include "xxx.ff/forcefield.itp"`
   改写成合并目录的**绝对路径**，所以合并目录要保留在磁盘上。

### 5.4 反向转换：GROMACS → AMBER（MM/PBSA 用）

gmx_MMPBSA 通过 **ParmEd** 把 GROMACS 拓扑转成 AMBER prmtop。
ParmEd 会用 `PBRadii` 指定的**连续介质半径集**重新赋予原子半径。

- **原生 AMBER 拓扑**：prmtop 里已存的 `RADII` / `SCREEN` 会被保留，`PBRadii` 不覆盖；
- **GROMACS 输入**：ParmEd 应用 `PBRadii` 生成 prmtop。
- 官方明确：*"原生 AMBER 与 GROMACS 转换两条路**不等价***，报同样的 advisory warning。

**这就是为什么推荐 AMBER 系蛋白 + GAFF 配体**：它离 ParmEd 的转换路径最短，
键合参数不容易出错（见第 0 节的实测提醒）。

---

## 6. gmx_MMPBSA 1.7.0 里与力场有关的参数

### 6.1 1.7.0 的默认值变化（相对 1.6 及更早）

| 变量 | 旧默认 | **1.7.0 默认** | 影响 |
|---|---|---|---|
| `igb` | 5 | **8** | 5 = GB-OBC2 → 8 = **GB-Neck2**（更新的模型）。旧脚本不写 `igb` 的话结果会变 |
| `PBRadii` | 3 | **4** | 3 = mbondi2 → 4 = **mbondi3**，与 `igb=8` 配套 |
| `pb exdi` | 80.0 | **78.5** | 溶剂介电常数 |

> **务必在方法学里写明 `igb` 和 `PBRadii` 的实际取值**，否则审稿人/同行无法复现。

### 6.2 `igb` 与 `PBRadii` 的配对（官方推荐表）

| `igb` | 模型 | 配套 `PBRadii` |
|---|---|---|
| 1 | GB-HCT | 2 (mbondi) |
| 2 | GB-OBC1 | 3 (mbondi2) |
| 5 | GB-OBC2 | 3 (mbondi2) |
| 7 | GB-Neck | 1 (bondi) |
| **8** | **GB-Neck2** | **4 (mbondi3)** ← 1.7.0 默认 |

`PBRadii` 也可以直接写名字：`PBRadii = mbondi3`。
官方特别注明：**`PBRadii` 只是"连续介质半径集"，它不选键合/非键力场**。

### 6.3 `source_force_field`：只是"来源标注"，不改任何东西

官方原文：

> `source_force_field` (Default = "auto") — Optional source force-field family override used for
> **provenance and advisories**. Allowed values are `auto`, `amber`, `charmm`, `opls`, `gromos`, `other`.
> This setting **never selects a radius set or changes a topology**.

白话：它只用来告诉程序"我这个体系是 AMBER 系的"，好让程序给出相应的提醒；
**它不会**替你切换半径集，**不会**改变拓扑。所以：

- 写 `source_force_field = "amber"` 是**正确且推荐**的（省得自动识别不出）；
- 但**不能**靠改它来"切换力场"——力场在 MD 拓扑里就已经定死了。

### 6.4 已废弃（legacy）的变量：不要写进新文件

官方原文：

> `forcefields`, `protein_forcefield`, `ligand_forcefield`, and `ions_parameters` — Legacy
> topology-preparation settings. They are retained only for compatibility with older input files and
> **are not used to parameterize current topology-based calculations. Do not add them to new input files.**

> Historical structure-to-topology workflows using `tleap`, `parmchk2`, and `-lm` are **not part of
> the current GROMACS calculation path**. A small-molecule ligand must already be represented in the
> supplied GROMACS topology.

**这两条要特别注意**：

1. 新建 `mmpbsa.in` **不要再写** `forcefields = "oldff/leaprc.ff99SB,leaprc.gaff"`。
2. 旧教程里的 `gmx_MMPBSA -lm ligand.mol2`（用 tleap 现算配体参数）**这条路已经不存在了**。
   现在必须**先**用 acpype/amb2gmx 把配体参数化好、合并进拓扑，再喂给 gmx_MMPBSA。

### 6.5 `idecomp`（能量分解方式）与 `print_res`

`idecomp` 取值（官方）：

| 值 | 模式 | 说明 |
|---|---|---|
| 1 | 逐残基 | 1-4 项计入内部能 |
| 2 | 逐残基 | 1-4 EEL 计入 EEL、1-4 VDW 计入 VDW |
| 3 | **成对** | 1-4 项计入内部能 |
| 4 | **成对** | 1-4 EEL/VDW 分别计入 EEL/VDW |

`print_res` 默认值是 **`"within 6"`**（1.7.0 起），意思是"打印受体和配体之间距离 6 Å 以内的残基"。

官方提醒：

- `print_res = "all"` —— *"This option is generally not recommended because most residue
  contributions are zero, so it can waste time and computational resources."*
  （本 skill 已从模板里改成 `within 6`。）
- `idecomp = 3/4`（成对模式）在大体系 + 大范围 print_res 时会产生**巨大的临时 mdout**
  （官方举例：250 个残基约 500 MB，因为是 62500 对），解析时内存占用也很大。

**其它格式要点**：

- `print_res` 支持 `within <距离>`、`"链/残基号"`（如 `"A/2-10,35,41 B/104"`）、`"all"`；
- 用 `within` 时若切出的残基 < 2 个，程序会自动把距离 +0.1 Å 直到 ≥ 2 个；
- **必须同时包含受体和配体的残基**，否则报
  `For decomposition analysis, you must define residues for both receptor and ligand!`；
- 带插入码的残基要单独写（如 `A/27B`），范围选择 `A/5-6B` 会报错。

### 6.6 其它与力场/溶剂相关的 1.7.0 变化

| 项 | 说明 |
|---|---|
| `qh_entropy` / `sander_apbs` | 1.7.0 **拒绝**执行，值必须保持 0。准谐振熵改用 NMODE / interaction entropy / C2 entropy |
| `explicit_waters` | 用 dASA mask 挑水；**TIP4P/OPC 这类多位点水在 sander 里会失败**，官方建议 `extra_points = "strip"`（并注明这会改变静电，适合受控相对比较），或改用三位点水（TIP3P/SPC） |
| `radii_audit = 1` | 输出逐原子半径来源 CSV；每次运行都会写 `GMXMMPBSA_radii.json` |
| `use_sander` | 1.7.0 中建议保持默认 1（用 sander 算 EEL） |
| `-cr` 参考结构 | 官方推荐用它保证残基在复合物/受体/配体三份拓扑间**完美对应** |

---

## 7. 推荐配方与完整命令

### 7.1 配方 A：蛋白-配体 + MM/PBSA（**本 skill 默认**）

```bash
# --- 蛋白 ---
gmx pdb2gmx -f receptor.pdb -o protein.gro -i posre.itp -p topol.top \
            -ff amber99sb -water tip3p -ignh < /dev/null

# --- 配体（GAFF）---
acpype -i ligand.mol2 -b ligand -a gaff -o gmx -c bcc -n 0

# --- 合并力场（GAFF 原子类型并进蛋白力场目录）---
python3 scripts/merge_gaff_ff.py \
  --lig-itp lig.acpype/lig_GMX.itp \
  --ff-src   /usr/local/share/gromacs/top/amber99sb.ff \
  --out-dir  ./amber99sb_gaff.ff \
  --strip-lig-itp ./lig.itp

# --- 组装复合物（拓扑里的力场引用会被改写成合并目录的绝对路径）---
python3 scripts/assemble_complex.py \
  --protein protein.gro --ligand lig.acpype/lig_GMX.gro \
  --topol topol.top --lig-itp ./lig.itp \
  --ff-dir ./amber99sb_gaff.ff --out-prefix complex

# --- 溶剂化（TIP3P 几何用 spc216.gro）---
gmx editconf -f complex.gro -o complex_c.gro -c -box 3.0 3.0 3.0
gmx solvate -cp complex_c.gro -cs /usr/local/share/gromacs/top/spc216.gro \
            -box 3.0 3.0 3.0 -o solv.gro -p topol.top
```

对应的 `mmpbsa.in`（**注意：没有 `forcefields` 那一行**）：

```fortran
&general
  sys_name            = "SystemName",
  startframe          = 1,
  endframe            = 50,
  interval            = 1,
  source_force_field = "amber",
  PBRadii             = 4,
  temperature         = 298.15,
  gmx_path            = "/usr/local/bin",
  solvated_trajectory = 0,
  verbose             = 2,
/

&gb
  igb         = 8,
  saltcon     = 0.150,
  intdiel     = 1.0,
  extdiel     = 78.5,
  surften     = 0.0072,
/

&decomp
  idecomp     = 2,
  print_res   = "within 6",
/
```

### 7.2 配方 B：CHARMM36

```bash
# 前提：已装 charmm36.ff（见 4.2）
gmx pdb2gmx -f receptor.pdb -o protein.gro -i posre.itp -p topol.top \
            -ff charmm36 -water tips3p -ignh < /dev/null
```

mdp 必须用 `assets/mdp/charmm36_*.mdp`，**不能**用通用模板。

MM/PBSA 侧：CHARMM 体系要用 `PBRadii = 7`（charmm_radii，见 6.2 之外的补充）

```fortran
&general
  source_force_field = "charmm",
  PBRadii             = 7,     ; charmm_radii
/
```

> `charmm_radii` 是 gmx_MMPBSA v1.5.0 加入的，官方注明：**只与用 CHARMM 力场制备的体系一起用**，
> 半径取自 CHARMM-GUI 的 PBEQ Solver 的 `radii.str`。

---

## 8. 常见错误速查表

| 现象 | 原因 | 处理 |
|---|---|---|
| `ERROR: Atomtype c3 not found` | GAFF 配体 atomtypes 没并进蛋白力场 | 用 `scripts/merge_gaff_ff.py` 合并 |
| `Fatal error: Syntax error in file topol.top` | 出现了两个 `[ atomtypes ]` 段 | 同上（合并而不是拼接） |
| `Could not find force field 'amber19sb'` | GROMACS 2025.4 没有这个目录 | 用 `amber99sb` 或 `amber99sb-ildn` |
| `Could not find water model 'opc'` | 未内置 OPC | 用 `tip3p`，或自行安装 OPC 文件 |
| `water model 'tips3p' not found` | 你在用 `amber99sb.ff` 而不是 `charmm27.ff`/`charmm36.ff` | `tips3p` 只在 CHARMM 目录里有 |
| CHARMM36 结果异常（密度/能量不合理） | mdp 没按官方改满 10 项 | 换用 `assets/mdp/charmm36_*.mdp` |
| grompp 报 `vdw-modifier` 相关警告 | 模板是通用模板（无 force-switch） | 同上 |
| `could not convert string to float: '*********'` | 某残基能量溢出（常见于 OPLS 拓扑转换） | 换 AMBER 系蛋白力场；见 `references/mmpbsa.md` |
| `IDECOMP cannot be 0` | 给了两个组即启用分解 | `idecomp = 1..4` |
| `For decomposition analysis, you must define residues for both receptor and ligand!` | `print_res` 只选了受体或只选了配体 | 用 `"within 6"`，或同时列出两边的残基 |
| gmx_MMPBSA 报 `forcefields` 相关警告 | 新文件里写了 legacy 变量 | 删掉 `forcefields` / `protein_forcefield` / `ligand_forcefield` / `ions_parameters` |
| MM/PBSA 结果与旧文章对不上 | `igb`/`PBRadii` 默认值在 1.7.0 变了（5→8、3→4） | 在方法学里显式写明 `igb` 与 `PBRadii` |
| 体系密度明显偏离文献值 | 用了 GROMOS + 单截断 | 官方已知问题，换 AMBER/CHARMM/OPLS |