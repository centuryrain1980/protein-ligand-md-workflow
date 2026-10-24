# 分子对接 (Vina / gnina)

蛋白-配体 MD 流程里，**配体必须先对接进口袋**，否则组装出来的复合物里配体是悬空的，
结合自由能没有物理意义。

## 目录

- [本机可用工具](#本机可用工具)
- [gnina 版本对照](#gnina-版本对照依据-github-官方-releases)
- [旧版 gnina (v1.1)](#旧版-gnina-v11)
- [gnina(推荐)](#gnina推荐)
- [Vina](#vina)
- [盒子中心怎么定](#盒子中心怎么定)
- [从对接结果到 MD 体系](#从对接结果到-md-体系)
- [踩过的坑](#踩过的坑)
- [本机实测结论](#本机实测结论)

---

## 本机可用工具

| 工具 | 版本 | 位置 | 状态 |
|---|---|---|---|
| **gnina** | **1.3.3 (master 6fe1ce2, 2026-06-30)** | `~/bin/gnina` → `~/apps/gnina/gnina`（包装脚本）→ `~/apps/gnina/gnina_v1.3.3` | ✅ 实测可对接（CPU） |
| gnina（旧版留档） | 1.1 (master e4cb380+, 2023-12) | `~/bin/gnina_v1.1` | ✅ 可运行，供结果对照 |
| **Vina** | 1.2.7 | `~/bin/vina` → `~/apps/vina/vina-bin` | ✅ 可运行 |
| meeko | 0.7+ | `~/.local/bin/mk_prepare_*.py` | ✅ 受体 PDBQT 制备 |
| Open Babel | 3.1.0 | `/usr/local/bin/obabel` | ✅ 格式转换 |

```bash
export PATH="$HOME/bin:$HOME/.local/bin:$PATH"
gnina --version        # gnina v1.3.3 master:6fe1ce2   Built Jun 30 2026.
gnina_v1.1 --version   # gnina v1.1 master:e4cb380+   Built Dec 18 2023.
vina --version         # AutoDock Vina v1.2.7
```

> **GPU**：gnina 支持 CUDA（`--device N`），但本机 WSL2 内**没有可用 GPU**
> （`/dev/dri`、`/dev/kfd` 都不存在；外接 AMD 显卡既没接入 WSL，CUDA 也不支持 AMD）。
> 跑 CPU 版即可（本机 22 线程，1231 原子受体、exhaustiveness=8 约 2 分钟）。

---

## gnina 版本对照（依据 GitHub 官方 Releases）

| 版本 | 发布日期 | 关键变化 | 官方提供的二进制 |
|---|---|---|---|
| **v1.3.3** | **2026-06-29** | 当前最新；在 v1.3 基础上的 minor bug fixes；初步实现 pygnina | `gnina.cuda12.8.static`（**仅此一个**） |
| v1.3.2 | 2025-07-08 | 提供 cuDNN 9 + CUDA 12.8 变体 | `gnina.1.3.2`（CPU）、`gnina.1.3.2.cuda12.8` |
| v1.3.1 | 2025-05-23 | bug 修复 | `gnina1.3.1` |
| **v1.3** | **2024-10-04** | **改用 PyTorch 实现 CNN；在 CrossDock2020 v1.3 上重训 CNN；引入知识蒸馏（knowledge distillation）打分函数** | `gnina`、`gnina.fix` |
| v1.1 | 2023-12-18 | 本 skill 之前的默认版本 | `gnina` |

### 升级建议

**结论：升到 v1.3.3。** 理由：

1. **CNN 模型更新**。v1.3 在 CrossDock2020 v1.3 上重训了 CNN 并用了知识蒸馏，
   对 pose 排序的精度优于 v1.1 的老模型；官方同时强调 v1.3 改用 PyTorch 实现，
   后续模型迭代更快。
2. **本机已实测跑通**（见下面的对比表），无回归问题。
3. 旧版 v1.1 保留在 `~/bin/gnina_v1.1`，用于复现旧结果或交叉验证。

### 三个必须注意的点

1. **不同版本的 CNN 分数不能横向比较**。CNN 模型换了（重训 + 蒸馏），
   `CNN pose score` / `CNN affinity` 的数值尺度和 v1.1 不是一回事，
   跨版本比大小没有意义；`affinity`（Vina 打分）不受影响。
   所以：**同一批配体要么全用新版本，要么全用旧版本。**
2. **官方对发布的二进制有免责声明**：Release notes 里写明发布的二进制
   *"for evaluation only, not production use"*（为兼容性牺牲了部分性能），
   正式生产建议自行从源码编译。本 skill 用的是官方二进制（CPU 模式足够本机用）。
3. **v1.3.3 只提供 CUDA 12.8 的静态构建**（`gnina.cuda12.8.static`，1.96 GB），
   即使跑 CPU 模式也**必须有 cuDNN 9 等 CUDA 运行时库**才能启动，
   见下面"CUDA 运行时依赖"。

### CUDA 运行时依赖（本机已配好）

v1.3.3 的静态构建在启动时仍要链接 `libcudnn.so.9` / `libcudart.so.12` /
`libcublas.so.12` / `libcusparse.so.12` 等。本机没有完整 CUDA Toolkit，
用 pip 装运行时库最省事（不需要管理员权限、不污染系统路径）：

```bash
pip3 install --break-system-packages --ignore-installed \
  -i https://pypi.tuna.tsinghua.edu.cn/simple \
  nvidia-cudnn-cu12 nvidia-cublas-cu12 nvidia-cuda-runtime-cu12 \
  nvidia-cusparse-cu12 nvidia-cusolver-cu12 nvidia-cufft-cu12 \
  nvidia-curand-cu12
```

装好后**必须**把库目录加进 `LD_LIBRARY_PATH`，否则报
`error while loading shared libraries: libcudnn.so.9: cannot open shared object file`。
本机用包装脚本 `~/apps/gnina/gnina` 自动处理（`~/bin/gnina` 指向它）：

```bash
#!/bin/bash
NV="$HOME/.local/lib/python3.14/site-packages/nvidia"
LIBS=$(find "$NV" -maxdepth 2 -type d -name lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="$LIBS${LD_LIBRARY_PATH:-}"
exec "$HOME/apps/gnina/gnina_v1.3.3" "$@"
```

> 直接在命令行 `export LD_LIBRARY_PATH=...` 也可以，但每开一个新终端都要重设，
> 所以统一走包装脚本。

### 本机实测对比（同一受体、同一配体、同一盒子、同一随机种子 42）

测试体系：`~/dock_test`（泛素受体 + 小分子配体），1231 原子受体，
盒子中心 (30.25, 28.91, 15.03)、边长 22 Å、exhaustiveness=4、num_modes=5、CPU 模式。

| 版本 | affinity (kcal/mol) | CNN pose score | CNN affinity |
|---|---|---|---|
| v1.1 | -3.66 / -2.28 / -3.41 / -2.55 / -2.40 | 0.687 ~ 0.833 | 2.87 ~ 3.09 |
| **v1.3.3** | **-2.97 / -4.30 / -4.10 / -3.28 / -4.78** | **0.771 ~ 0.875** | 3.15 ~ 4.38 |

耗时（另一组 exhaustiveness=8、同一台机器、连续两次实测）：

| 版本 | 耗时 |
|---|---|
| v1.1 | 122 秒 |
| **v1.3.3** | **111 秒**（与 v1.1 同量级，略快） |

**怎么读这张表**：

- **affinity 整体更负**（v1.3.3 最优 -4.78 vs v1.1 最优 -3.66）——新版找到了更好的构象。
  **但这不等于"新版比旧版准"**：同一种子的两次搜索本身就可能不同，
  严格比较需要多个不同种子的重复实验。
- **CNN pose score 整体更高**（0.771 ~ 0.875 vs 0.687 ~ 0.833）——新版 CNN 对这些 pose 评价更高。
- **CNN affinity 的绝对值完全不同**（3.1 ~ 4.4 vs 2.9 ~ 3.1）——这正是上面第 1 条提醒的证据：
  **换了模型，分数尺度就变了，跨版本不能比**。
- 有 NVIDIA GPU 时加 `--device 0` 会快一个数量级。

---

## 旧版 gnina (v1.1)

保留在 `~/bin/gnina_v1.1`，用于复现旧结果。命令行参数与 v1.3.3 完全兼容：

```bash
gnina_v1.1 --receptor receptor.pdb --ligand ligand.pdbqt \
      --center_x 30.25 --center_y 28.91 --center_z 15.03 \
      --size_x 22 --size_y 22 --size_z 22 \
      --exhaustiveness 4 --num_modes 5 --no_gpu \
      --out docked.sdf --log gnina.log
```

---

## gnina(推荐)

gnina = Vina 打分 + **CNN 深度学习重打分**，对小分子结合 pose 的排序通常比纯 Vina 更准。

**关键优势：直接吃 PDB 受体**，内部自动转 PDBQT，不需要 meeko。

```bash
gnina --receptor receptor.pdb \\
      --ligand   ligand.pdbqt \\
      --center_x $CX --center_y $CY --center_z $CZ \\
      --size_x 22 --size_y 22 --size_z 22 \\
      --exhaustiveness 32 --num_modes 9 \\
      --no_gpu \\
      --out docked.sdf --log gnina.log
```

输出 `docked.sdf` 含多个构象（用 `$$$$` 分隔），日志末尾给出：

```
mode |  affinity  |  intramol  |    CNN     |   CNN
     | (kcal/mol) | (kcal/mol) | pose score | affinity
-----+------------+------------+------------+----------
    1       -5.83       -0.42       0.7717      2.984
```

- `affinity`：Vina 的对接打分（kcal/mol，越负越好）
- `CNN pose score` / `CNN affinity`：CNN 重打分，**挑 pose 时优先看这个**
- `intramol`：配体自身内部能量，异常正值说明构象扭曲

常用参数：

| 参数 | 说明 |
|---|---|
| `--autobox <ligand>` | 自动以配体为中心定盒子（盲对接用） |
| `--cnn_scoring rescore` | 只用 CNN 重打分已有构象 |
| `--minimize` | 对输出构象做能量最小化 |
| `--seed <n>` | 随机种子 |
| `--local_only` | 只在配体附近微调 |
| `--cspg` / `--cgnn` | 自定义 CNN 模型 |

---

## Vina

**必须有 PDBQT 受体**（Vina 不接受 PDB 受体）。

```bash
# 1) 受体 PDBQT (meeko)
mk_prepare_receptor.py -i receptor.pdb -o receptor.pdbqt -p 7.4

# 2) 对接
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \\
     --center_x $CX --center_y $CY --center_z $CZ \\
     --size_x 22 --size_y 22 --size_z 22 \\
     --exhaustiveness 32 --num_modes 9 \\
     --out docked.pdbqt --log vina.log
```

配体 PDBQT（Open Babel 3.x：`-o` 接格式，`-O` 接文件名）：

```bash
obabel -imol2 -opdbqt -Oligand.pdbqt ligand.mol2
```

> **Open Babel 3.x 转不了蛋白成 PDBQT**（`0 molecules converted`）。
> 受体 PDBQT 用 meeko；配体 PDBQT 用 obabel 或 meeko 都可以。
>
> **实测提醒**：本机 meeko 对 GROMACS 导出的 PDB 会一直报
> `Use --wanted_altloc to set variants for specific residues.`（尽管 PDB 里没有 altLoc 记录），
> prody 与 gemmi 两个读取器都一样。**因此 Vina 路径在本机不如 gnina 顺**——
> gnina 直接接受 PDB 受体并内部完成同样的制备。
>
> 需要 Vina 的 PDBQT 受体时，可行方案：
> 1. 用 AutoDockTools 的 `prepare_receptor4.py`（需 conda 装 `autodock-vina` 配套环境）
> 2. 用 `pdb2pqr` 生成 PQR 后再转
> 3. 换用结构更规整的受体 PDB（如去掉 TER/ANISOU 记录）

---

## 盒子中心怎么定

**这是对接结果是否合理的关键**。盒子没对准口袋时，gnina/Vina 依然会输出构象，
但亲和力全是 0.00（配体落在蛋白外面，没有任何相互作用）——这很容易被误读成"结合很弱"。

三种定法：

```bash
# 1) 已知口袋: 用结构信息里的口袋残基中心
#    (从 PDB 的 SITE/HETATM 记录, 或从对接软件的输出)

# 2) 已知共晶配体: 用它的位置
obabel -ipdb -oxyz -Olig_center.xyz cocrystal_ligand.pdb

# 3) 盲对接: 用蛋白中心 + 大盒子
python3 - <<'PY'
L=open("receptor.pdb").read().splitlines()
c=[(float(l[30:38]),float(l[38:46]),float(l[46:54])) for l in L if l.startswith("ATOM")]
n=len(c)
print(sum(x for x,_,_ in c)/n, sum(y for _,y,_ in c)/n, sum(z for _,_,z in c)/n)
PY
```

盒子边长经验值：口袋直径 + 8~10 Å，即 20~25 Å。

---

## 从对接结果到 MD 体系

对接输出的是 PDBQT/SDF，**不能直接进 GROMACS**，需要：

1. 取最优构象 → SDF/PDB
2. 重新加氢、确定质子化态（`obabel`/pdb2pqr/PropKa）
3. 用 acpype 生成 GAFF 参数
4. **把配体坐标并入蛋白坐标**（`scripts/assemble_complex.py` 负责，含残基号重编）
5. 溶剂化、加离子、平衡、生产 MD

`scripts/setup_ligand.sh` 覆盖了第 1–4 步（SMILES 起手）。

---

## 踩过的坑

| 现象 | 原因 | 处理 |
|---|---|---|
| `error while loading shared libraries: libcudnn.so.9` | 直接运行 v1.3.3 原始二进制，没带 CUDA 运行时库路径 | 用包装脚本 `~/apps/gnina/gnina`（`~/bin/gnina`），它会自动设 `LD_LIBRARY_PATH` |
| `affinity` 全是 **0.00** | 盒子没对准蛋白（或受体 PDBQT 是空的） | 核对盒子中心；检查受体 PDBQT 原子数 > 0 |
| Vina 打印 usage 摘要就退出 | 受体不是 PDBQT | 用 meeko 制备 |
| obabel 转蛋白 PDBQT 得到 **0 原子** | Open Babel 3.x 不支持蛋白→PDBQT | 用 meeko |
| meeko 报 `Element '' not found` | 输入 PDB 没有元素列 | 用 `gmx editconf`/`trjconv` 导出的规范 PDB，别手工裁列 |
| `mk_prepare_receptor.py` 报缺 prody/gemmi/rdkit | meeko 依赖链 | `pip install meeko rdkit gemmi prody` |
| meeko 报 `Use --wanted_altloc to set variants` 但 PDB 里没有 altLoc | meeko/prody 对 GROMACS 写的 PDB 的检测问题 | **直接用 gnina**（吃 PDB 受体，无需 PDBQT）；或用 AutoDockTools 的 `prepare_receptor4.py` 产出 PDBQT |
| 配体 PDBQT 原子数比配体少 | PDBQT 默认丢非极性氢 | 正常；对接时由 Vina 自动补 |
| `~/bin/vina` 与 `~/bin/gnina` 不在 PATH | `~/bin` 只在登录 shell 加入 | `export PATH="$HOME/bin:$PATH"` |
| 换 gnina 版本后 CNN 分数与旧文章对不上 | v1.3 重训并蒸馏了 CNN 模型，数值尺度变了 | 不要跨版本比较；同一批配体统一版本 |

---

## 本机实测结论

| 路径 | 状态 | 说明 |
|---|---|---|
| **gnina 1.3.3 + PDB 受体 + PDBQT 配体** | ✅ 完全跑通 | 1231 原子受体、盒对准蛋白中心 → 5 个构象、affinity -2.97 ~ -4.78 kcal/mol、CNN pose score 0.7714 ~ 0.8749 |
| gnina 1.1（旧版留档） | ✅ 可运行 | 同一体系 affinity -3.66 ~ -2.28、CNN pose score 0.6870 ~ 0.8327 |
| Vina + PDB 受体 | ❌ | Vina 打印 usage 摘要退出，必须给 PDBQT 受体 |
| Vina + meeko PDBQT 受体 | ⚠️ | meeko 在本机对该 PDB 报 altLoc 错误，见上 |

**结论：优先用 gnina。** 它既是 Vina 的打分 + CNN 重打分，又省掉了受体 PDBQT 制备这一环。
默认用 `~/bin/gnina`（v1.3.3），需要复现旧结果时用 `~/bin/gnina_v1.1`。

一键脚本：`bash scripts/dock_ligand.sh <受体PDB> <配体mol2> [cx cy cz] [size]`
