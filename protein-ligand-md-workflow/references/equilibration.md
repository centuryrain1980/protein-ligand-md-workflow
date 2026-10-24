# 平衡与生产 MD

## 目录

- [阶段设计](#阶段设计)
- [1. 能量最小化 EM](#1-能量最小化-em)
- [2. NVT 平衡](#2-nvt-平衡)
- [3. NPT 平衡](#3-npt-平衡)
- [4. 生产 MD](#4-生产-md)
- [5. 位置限制的逐步放开](#5-位置限制的逐步放开)
- [6. 性能调优](#6-性能调优)
- [7. 排错](#7-排错)

---

## 阶段设计

标准四阶段协议，每阶段解决一类问题：

| 阶段 | 解决什么 | 典型时长 |
|---|---|---|
| EM | 消除加离子/溶剂化带来的严重接触 | 到收敛即可（本机 14346 步） |
| NVT | 升温到目标温度，松开部分约束 | 100 ps – 1 ns |
| NPT | 放开体积，压到真实密度 | 100 ps – 1 ns |
| Production | 无约束的长时间采样 | 10 ns – 数百 ns |

**约束的逐步放开**（`define = -DPOSRES`）：

1. EM：无位置限制
2. NVT：位置限制蛋白（`define = -DPOSRES`）
3. NPT：位置限制蛋白+配体（`define = -DPOSRES_FC` 或自定义）
4. Production：无位置限制（去掉 `define`）

这样能避免在剧烈涨落过程中把蛋白结构扯坏。

---

## 1. 能量最小化 EM

用 steepest descent 把体系带到局部能量极小。

```bash
gmx grompp -f mdp/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -s em.tpr -deffnm em -ntmpi 1 -ntomp 6 -pin on
```

参数：

```
integrator     = steep
emtol          = 10.0      ; kJ/mol/nm，力的收敛阈值
emstep         = 0.01
nsteps         = 50000
nstlist        = 10
cutoff-scheme  = Verlet
coulombtype    = PME
rcoulomb       = 1.0
rvdw           = 1.0
pbc            = xyz
```

**收敛判据**（md.log）：

```
Steepest Descents converged to Fmax < 10 in N steps
Steepest Descents converged to machine precision in N steps   <- 更好
```

没收敛（跑满 50000 步）通常说明体系里有严重冲突。检查：
加离子后有没有重新最小化、溶剂化时有没有异常重叠、盒子里是不是漏了配体。

---

## 2. NVT 平衡

恒温恒容，把体系从初始构象松弛到目标温度。

```bash
gmx grompp -f mdp/nvt.mdp -c em.gro -r solv_ions.gro -p topol.top -o nvt.tpr
gmx mdrun -s nvt.tpr -deffnm nvt -ntmpi 1 -ntomp 6 -pin on
```

**关键点**

- `-r solv_ions.gro` 必须传！mdp 里有 `define = -DPOSRES`，
  缺 `-r` 报 `Cannot find position restraint file restraint.gro (option -r)`。
- `gen_vel = yes` + `gen_temp = 300` 生成初始速度。
- `tcoupl = V-rescale`（修正 Berendsen 恒温器），比老式 Berendsen 更适合平衡与生产。
- `tc-grps = Water Non-Water` 分组耦合：水的比热大，与溶质分开耦合更准确。
- `tau_t = 0.1` ps 是时间常数，0.1 ps 属"中等耦合"。

**2025 版必须改的设置**：

```
comm-mode               = Linear      ; 不是 Angular
comm-grps               = System      ; 不是 Protein
```

否则 grompp 会因 2 个警告而失败（详见 `gromacs-2025-api.md`）。

---

## 3. NPT 平衡

恒压恒温，让盒子涨缩到真实密度。

```bash
gmx grompp -f mdp/npt.mdp -c nvt.gro -r solv_ions.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -s npt.tpr -deffnm npt -ntmpi 1 -ntomp 6 -pin on
```

| 选项 | 值 | 说明 |
|---|---|---|
| `pcoupl` | `Parrinello-Rahman` | 平衡阶段的标准选择 |
| `pcoupltype` | `isotropic` | 各向同性；膜体系用 `semiisotropic` |
| `tau_p` | `2.0` ps | 压力耦合时间常数 |
| `ref_p` | `1.0` bar | 参考压力 |
| `compressibility` | `4.5e-5` bar⁻¹ | 水的等温压缩系数 |
| `refcoord_scaling` | `com` | 随盒子一起缩放 |
| `continuation` | `yes` | 延续 NVT 的速度与检查点 |
| `gen_vel` | `no` | 不重新生成速度 |

**必须传 `-t nvt.cpt`**（延续 NVT 检查点），否则 NVT/NPT 之间有能量跳变。

**收敛判据**：md.log 里密度与压力应趋于平稳，密度收敛到约 0.98–1.05 g/cm³（取决于力场）。

---

## 4. 生产 MD

```bash
gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -s md.tpr -deffnm md -ntmpi 1 -ntomp 6 -pin on -x md.xtc
```

生产阶段的 mdp 相对 NPT 的改动：

- 去掉 `define = -DPOSRES`（不再限制位置）
- `nstxout-compressed = 5000`（每 10 ps 存一帧压缩轨迹，xtc 体积小）
- `nstxout = 0`、`nstvout = 0`（不存全精度轨迹和速度，节省磁盘）
- `nstenergy = 5000`
- **加能量组**（若要做结合能分析）：
  ```
  energygrps             = Protein LIG Protein_LIG
  ```
  这会让 .edr 里出现 `Coul-SR:Protein-LIG` 之类的组间能量项。

**采样频率权衡**：

| 用途 | 建议 `nstxout-compressed` |
|---|---|
| 结合自由能 MM/PBSA | 100–200（每 0.2–0.4 ps） |
| 常规构象分析 | 5000–10000（每 10–20 ps） |
| 氢键占有率 | 100–500 |
| 自由能面 (FEL) | 100–1000 |

**多重复模拟**：单条轨迹可能陷入局部极小。做 3–5 条独立重复（不同 `gen_seed`），
先各自平衡再合并成concatenated 轨迹，用于充分采样。

---

## 5. 位置限制的逐步放开

如果担心放开约束后结构崩掉，用分阶段放开：

```
; NVT 阶段
define                  = -DPOSRES        ; 蛋白强限制
define                  = -DPOSRES_FC     ; 柔性约束版本

; 生产阶段
; （不写 define，或注释掉）
```

更精细的做法是自定义位置限制：只限制蛋白的 Cα，或只限制结合位点残基。
用 `gmx make_ndx` 建组，再在 topol.top 里加：

```
[ position_restraints ]
; atomfunct  ai  aj  type  fc  fcstep
1           1   2   1     1000 1000
```

---

## 6. 性能调优

### 本机（Intel Core Ultra 9 185H，11 物理核 / 22 线程，混合架构）

```bash
# 单体系（推荐）：只用 6 个大核，绑定线程
gmx mdrun -s md.tpr -deffnm md -ntmpi 1 -ntomp 6 -pin on

# 大体系：MPI 多进程
mpirun -np 4 gmx_mpi mdrun -s md.tpr -deffnm md
```

- `-ntomp 6` 配合 `-pin on` 比用满 22 线程快（能效核拖慢，且线程争抢）。
- GROMACS 会警告 "The number of threads is not equal to the number of (logical) cpus and
  the -pin option is set to auto" —— 显式加 `-pin on` 消除。
- 实测参考：2731 原子水盒，4 线程 ~248 ns/day（EM/NVT），6 线程 ~221–228 ns/day。

### 体系规模与线程数

| 原子数 | 建议 |
|---|---|
| < 1 万 | `-ntomp 4 ~ 6` |
| 1–10 万 | `-ntomp 8 ~ 11` |
| > 10 万 | MPI：`mpirun -np 4 gmx_mpi mdrun`（注意 `gmx_mpi` 不接受 `-ntmpi`） |

### 检查点与断点续跑

```bash
# 从检查点继续
gmx mdrun -s md.tpr -c md.gro -t md.cpt -deffnm md -append
```

---

## 7. 排错

| 现象 | 原因 | 处理 |
|---|---|---|
| `Too many warnings (N)` | 2025 版警告变错误 | 见 `gromacs-2025-api.md` 第四节 |
| `Cannot find position restraint file` | 有 `-DPOSRES` 但没传 `-r` | grompp 加 `-r <坐标>` |
| `Group XXX referenced in the .mdp file was not found` | 传 `-n` 后默认组被替换 | 自定义 ndx 里补上该组，或不传 `-n` |
| `N atoms are not part of any of the T-Coupling groups` | TC 组没覆盖全部原子 | 用 make_ndx 的默认 `Water`/`non-Water` |
| `number of coordinates ... does not match topology` | `[ molecules ]` 数量对不上 | 用 `solvate -p`/`genion -p` 自动更新，或手工核对 |
| 没有生成 .xtc | mdp 缺 `nstxout-compressed` | 补上，并给 mdrun 加 `-x md.xtc` |
| 温度忽高忽低 | 耦合太强 / 时间常数太小 | `tau_t` 调大到 0.5–1.0 |
| 密度持续漂移 | NPT 未收敛 | 延长 NPT，或检查 `refcoord_scaling` |
| 蛋白结构在平衡期就散架 | 位置限制放开太早 | 延长带约束的 NVT，或减小 `tau_t` |
| mdrun 跑得极慢 | 线程数超配 | 降到 `-ntomp 6` 并加 `-pin on` |
