# GROMACS 2025 命令变更与排错

本文所有条目都在本机 GROMACS 2025.4 上**实测确认**。旧教程（2023/2024）与本版差异很大，直接照抄会失败。

---

## 目录

- [一、会静默卡死的命令（最危险）](#一会静默卡死的命令最危险)
- [二、参数语义变更](#二参数语义变更)
- [三、已移除的命令与选项](#三已移除的命令与选项)
- [四、警告升级为错误](#四警告升级为错误)
- [五、交互式选择的标准喂法](#五交互式选择的标准喂法)
- [六、索引组的规则](#六索引组的规则)
- [七、输出文件命名的坑](#七输出文件命名的坑)
- [八、快速对照表](#八快速对照表)

---

## 一、会静默卡死的命令（最危险）

### `gmx pdb2gmx` —— 缺 `-ff` 或 `-water` 会无限循环

GROMACS 2025 的 pdb2gmx 会依次要求选择**力场**和**水模型**。在非交互环境（无 TTY、stdin 为空）下，
它读完 EOF 后进入用户态死循环：**占用 100% CPU、不打印任何错误、不会退出**。

```bash
# 错误：会永久卡住（无任何报错）
gmx pdb2gmx -f protein.pdb -o out.gro

# 正确：显式给出力场和水模型，全程无交互
gmx pdb2gmx -f protein.pdb -o protein.gro -i posre.itp -p topol.top \
            -ff oplsaa -water spc -ignh < /dev/null
```

- `-ff` 可选值：`oplsaa`、`amber99sb-ildn`、`charmm36-jul2021` 等
- `-water` 可选值：`spc`、`tip3p`、`tip4p`、`tip5p`、`none`
- `-ignh` 忽略输入中的氢原子（PDB 通常没有，无需担心）
- `-ter` 会额外询问末端质子化状态，非交互环境**不要**加

> 诊断手法：`ps` 看 CPU 占用 100%，再 `cat /proc/<pid>/syscall` 与 `ls -l /proc/<pid>/fd`，
> 若看到 `/dev/ptmx` 且 `read(0)` 返回 0，就是这类死循环。

### `gmx editconf`、`gmx rms`、`gmx trjconv` 等带原子组选择的工具

这些工具会提示选择原子组。无输入时行为不一致（有的报错、有的循环）。
**统一做法**：用 `printf '组名\n组名\n' | gmx ...` 显式喂选择，或 `< /dev/null` 接受默认值。

---

## 二、参数语义变更

### `gmx editconf` 的 `-center` 从布尔开关变成了向量

老写法 `gmx editconf -f x.gro -o y.gro -center -box 3 3 3` 在 2025 版报
`In command-line option -center: Too few (valid) values`。

2025 版语义：

| 选项 | 2025 语义 |
|---|---|
| `-[no]c` | 布尔：把分子居中到盒子里（由 -box / -d 隐含） |
| `-center <vector>` | 三维向量：把几何中心平移到 (x,y,z) |
| `-align <vector>` | 三维向量 |
| `-scale <vector>` | 三维向量（取代旧的 `-d <real>` 标量缩放） |
| `-box <vector>` | 三维向量 |

```bash
# 2025 版：居中 + 建盒（-c 是布尔）
gmx editconf -f protein.gro -o protein_c.gro -c -box 3.0 3.0 3.0
```

editconf 也不再接受原子组参数（旧的 `-c <group>` 语义已变）。

### `gmx grompp` 不再自动选择力场

`-ff`（力场目录）、`-water`（水模型）、`-i`（分子信息）三个选项**全部被移除**。
必须自己写 `topol.top`，其中 `#include` 力场文件，并显式声明 `[ molecules ]` 段的分子数量。

缺 `[ molecules ]` 数量时报：`No molecules were defined in the system`。

### `gmx energy` 选的是能量项编号，不是原子组

老教程里 `gmx energy` 先问原子组，2025 版**只问能量项**，用编号或名称选择，空行结束。

完整的能量项列表（63 项，本机实测）：

```
 1 Bond          2 Angle         3 Proper-Dih.    4 Ryckaert-Bell.
 5 LJ-14         6 Coulomb-14    7 LJ-(SR)        8 Disper.-corr.
 9 Coulomb-(SR) 10 Coul.-recip. 11 Potential     12 Kinetic-En.
13 Total-Energy 14 Conserved-En. 15 Temperature  16 Pres.-DC
17 Pressure    18 Constr.-rmsd  19 Box-X         20 Box-Y
21 Box-Z        22 Volume       23 Density       24 pV
25 Enthalpy    26 Vir-XX       27 Vir-XY        28 Vir-XZ
29 Vir-YX      30 Vir-YY       31 Vir-YZ        32 Vir-ZX
33 Vir-ZY      34 Vir-ZZ       35 Pres-XX       36 Pres-XY
37 Pres-XZ     38 Pres-YX      39 Pres-YY       40 Pres-YZ
41 Pres-ZX     42 Pres-ZY      43 Pres-ZZ       44 #Surf*SurfTen
45 Box-Vel-XX  46 Box-Vel-YY   47 Box-Vel-ZZ
48 Coul-SR:Protein-Protein      49 LJ-SR:Protein-Protein
50 Coul-14:Protein-Protein      51 LJ-14:Protein-Protein
52 Coul-SR:Protein-rest         53 LJ-SR:Protein-rest
54 Coul-14:Protein-rest         55 LJ-14:Protein-rest
56 Coul-SR:rest-rest            57 LJ-SR:rest-rest
58 Coul-14:rest-rest            59 LJ-14:rest-rest
60 T-Water      61 T-non-Water  62 Lamb-Water    63 Lamb-non-Water
```

其中 **48–59 是组间能量**，只要 mdp 里设了 `energygrps`，就能直接拿到蛋白-配体相互作用能，
比"整体相减"精确得多（没有大数相消误差）。

`-g` 选项（输出两组的相互作用能）已被移除。

### `gmx solvate` / `gmx genion` 选项

`solvate` 用 `-scale`（标量，缩放溶剂盒）而不是别的；`genion` 用 `-pname/-nname` 指定离子名。

---

## 二点五、更多 2025 版行为变更

### `gmx trjconv`：PBC 处理与旋转拟合必须分两步

同一条命令里同时给 `-pbc` 和 `-fit` 会直接报错：

```
PBC condition treatment does not work together with rotational fit.
Please do the PBC condition treatment first and then run trjconv in a
second step for the rotational fit.
```

正确写法（两步）：

```bash
# 第一步: 处理周期性边界
printf 'Protein\nProtein\n' | gmx trjconv -s md.tpr -f md.xtc \
    -o whole.xtc -pbc whole

# 第二步: 居中 + 旋转拟合
printf 'Protein\nProtein\nProtein\n' | gmx trjconv -s md.tpr -f whole.xtc \
    -o traj_fit.xtc -center -fit rot+trans
```

**注意 trjconv 的组选择次数**（喂少了会读到 EOF 而失败）：

| 命令 | 需要的选择次数 |
|---|---|
| `-pbc whole` | 2（PBC 组 + 输出组） |
| `-center -fit rot+trans` | 3（拟合组 + 居中组 + 输出组） |

### `gmx hbond` 被重写

2025 版 `gmx hbond` 的选项语义全变：

| 选项 | 旧版 | 2025 版 |
|---|---|---|
| `-h` | 氢键图谱的输出文件名 | **布尔开关**（合并氢键对），传文件名报 `Invalid value` |
| `-o` | 输出 ndx | 输出 `.ndx` |
| `-hbm` | 无 | 氢键图谱 `.xpm`（在 `hbond-legacy` 里） |
| `-hbr` | 无 | 氢键距离截断，默认 0.35 nm |
| `-fgroup/-r/-t` | 无 | 新的 selection 语法 |

```bash
# 新版: 统计
printf 'Protein\nProtein\n' | gmx hbond -s md.tpr -f traj.xtc \
    -num hbond_num.xvg -dist hbond_dist.xvg

# 旧版: 氢键存在性图谱 (xpm)
printf 'Protein\nProtein\n' | gmx hbond-legacy -s md.tpr -f traj.xtc \
    -num hbond_num2.xvg -hbm hbond.xpm
```

`gmx hbond-legacy` 在 2025 版仍然提供，用于生成旧式的氢键矩阵 xpm。


### `gmx editconf` 在索引含越界组时会段错误

若 `-n` 指定的索引文件里有原子序号超出结构原子数的组，GROMACS 2025 的 `editconf`
**直接 Segmentation fault**（退出码 139），不打印任何错误信息。

常见触发场景：给一个 21 原子的配体结构配了一份含 6541 原子组的索引。

```bash
# 错误做法: 索引写复合物编号, 结构是配体
gmx editconf -f lig.pdb -n complex_index.ndx -o out.pdb    # -> 段错误

# 正确: 索引用该结构自身的局部编号
gmx editconf -f lig.pdb -n lig_local.ndx -o out.pdb
```

诊断：`ls -l /proc/<pid>/fd` 只能看到文件，要用 gdb/strace 才能看到崩溃点；
更简单的预防办法是**每次都检查索引组的最大原子号是否 <= 结构的原子数**。

## 三、已移除的命令与选项

| 移除项 | 替代方案 |
|---|---|
| `gmx genbox` | `gmx solvate -cp 溶质.gro -cs spc216.gro -box X X X -o out.gro -p topol.top` |
| `gmx tprconv` | 无直接替代；用 `gmx trjconv` |
| `gmx do_dssp` | `gmx dssp`（2023+ 输出符号矩阵，用 `dit dssp` 转换） |
| `grompp -ff` / `-water` / `-i` | 自写 `topol.top` + `[ molecules ]` 段 |
| `gmx energy -g` | `energygrps` + 选 48–59 号组间能量项 |

`solvate` 会**追加**而非替换 `topol.top` 里的 `[ molecules ]`：

```
[ molecules ]
; Compound        #mols
Protein_chain_A     1
SOL                58      <- pdb2gmx 写入的
SOL               442      <- solvate 追加的
```

这是合法的，grompp 会正确处理。但如果你手动改过这个段，注意别把行删错。

---

## 四、警告升级为错误

GROMACS 2025 起，**多个警告会直接终止 grompp**：

```
Fatal error:
Too many warnings (2).
```

最常见的两个警告来自中心-of-mass 设置：

```
WARNING 1: Removing the rotation around the center of mass in a periodic system,
            this can lead to artifacts.
WARNING 2: 1500 atoms are not part of any center of mass motion removal group.
            In most cases one should use one group for the whole system.
```

**正确修法**（不要用 `-maxwarn` 掩盖）：

```
comm-mode               = Linear
comm-grps               = System
```

周期体系里整体旋转的伪影问题，正确做法是不要用 `Angular` 模式移除旋转，
或只做平移修正（`Linear`）。把 Protein/配体等分组的 COM 移除用于周期体系是不安全的。

真正无害的警告才用 `-maxwarn 1` 放行。

---

## 五、交互式选择的标准喂法

| 场景 | 命令 |
|---|---|
| `gmx rms` / `rmsf` / `gyrate` / `sasa` / `hbond` | `printf 'Protein\nProtein\n' \| gmx rms -s t.tpr -f traj.xtc -o rms.xvg` |
| `gmx trjconv` | `printf 'System\nSystem\n' \| gmx trjconv -s t.tpr -f traj.xtc -o out.xtc` |
| `gmx genion` | `printf 'SOL\n' \| gmx genion -s ions.tpr -o solv.gro -p topol.top -pname NA -nname CL -neutral` |
| `gmx energy` | `printf '7\n8\n9\n10\n' \| gmx energy -f md.edr -o pro.xvg`（选能量项编号） |
| `gmx make_ndx` | `printf 'q\n' \| gmx make_ndx -f conf.gro -o index.ndx`（只生成默认组） |
| `gmx pdb2gmx` | 给 `-ff` 和 `-water`，不靠交互 |
| `gmx mdrun` | `gmx mdrun -s t.tpr -deffnm md`（无交互） |

`mdrun` 的 `-o` 选项在 2025 版**只接受 `.trr`/`.cpt`/`.tng`**，要 `.gro` 就用 `-deffnm`（会自动写 `<name>.gro`）。

---

## 六、索引组的规则

### 默认组名（`gmx make_ndx` 自动生成）

蛋白-配体 + 溶剂体系下典型默认组：

```
System, Protein, Protein-H, C-alpha, Backbone, MainChain, MainChain+Cb,
MainChain+H, SideChain, SideChain-H, Prot-Masses, non-Protein, Water,
Water_and_ions, SOL, non-Water, Ion, ...
```

注意是 **`non-Water`（小写 n）**。

### 组名匹配不区分大小写

mdp 里写 `tc-grps = Water Non-Water` 能匹配到默认组 `non-Water`，已实测确认。

### 传了 `-n` 就完全替换默认组

```bash
gmx grompp -f mdp.mdp -c conf.gro -p topol.top -n index.ndx -o out.tpr
```

用了 `-n` 之后，**默认组全部消失**，自定义文件里必须包含 mdp 引用的每一个组名，
否则报 `Group XXX referenced in the .mdp file was not found in the list of index groups`。

### TC 耦合组必须覆盖全部原子

`tc-grps = Water Non-Water` 要求两个组**不重不漏地覆盖全部原子**，
否则报 `N atoms are not part of any of the T-Coupling groups`。

用 `gmx make_ndx` 生成默认组时，`Water` 与 `non-Water` 天然互补，最省事。

---

## 七、输出文件命名的坑

### `grompp` 的输出名取拓扑文件名前 5 个字符

```
topol.top        -> topol.tpr
topol_big.top    -> topol.tpr     <-- 撞名！
```

**不同体系务必显式指定 `-o`**：

```bash
gmx grompp -f mdp.mdp -c conf.gro -p topol.top -o big.tpr
```

### 坐标/拓扑文件格式规范

`gmx pdb2gmx` 对输入 PDB 的原子命名很严格。残基缺主链原子时报：

```
Residue 1 named ALA of a molecule in the input file was mapped
to an entry in the topology database, but the atom C used in
that entry is not found in the input file.
```

---

## 八、快速对照表

| 操作 | 旧写法（≤2024） | GROMACS 2025 正确写法 |
|---|---|---|
| pdb2gmx 选力场 | 交互选择 | `-ff oplsaa -water spc` |
| 居中建盒 | `editconf -center -box 3 3 3` | `editconf -c -box 3.0 3.0 3.0` |
| 指定力场到 grompp | `grompp -ff oplsaa -water spc` | 自写 topol.top + `#include` |
| 声明分子数 | `grompp -i coord.gro` | topol.top 里写 `[ molecules ]` |
| 输出 tpr | `grompp -o out.tpr` | 同左（必须 .tpr 后缀） |
| 填水 | `genbox -f x.gro -o y.gro -box 5 5 5` | `solvate -cp x.gro -cs spc216.gro -box 5 5 5 -o y.gro -p topol.top` |
| 压缩轨迹 | `nstxout = 500` | `nstxout-compressed = 500` + `-x traj.xtc` |
| 能量选择 | 选原子组 | 选能量项编号（如 `7,8,9,10`） |
| 二级结构 | `do_dssp` | `gmx dssp` + `dit dssp` |
| COM 移除 | `comm-mode = Angular` + `comm-grps = Protein` | `comm-mode = Linear` + `comm-grps = System` |
| 警告处理 | 忽略 | 消除警告，或 `-maxwarn N` |
| `mdrun` 输出 gro | `mdrun -o out.gro` | `mdrun -deffnm md`（自动写 md.gro） |
