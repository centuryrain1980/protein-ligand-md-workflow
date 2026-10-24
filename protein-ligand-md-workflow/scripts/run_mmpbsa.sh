#!/bin/bash
# 一键运行 gmx_MMPBSA 结合自由能计算(MM/GBSA)
# 用法: bash run_mmpbsa.sh <工作目录> [受体组] [配体组] [起始帧] [结束帧]
set -u
export MPLBACKEND=Agg

WORK="${1:?用法: run_mmpbsa.sh <工作目录> [受体组] [配体组] [startframe] [endframe]}"
PROT="${2:-Protein}"
LIG="${3:-LIG}"
FSTART="${4:-1}"
FEND="${5:-50}"

# ---- 力场相关参数(可按体系覆盖) ----
# SRC_FF : source_force_field, 只用于来源标注与提示, 不会改变拓扑或半径集
#          取值 auto|amber|charmm|opls|gromos|other   官方文档: 见 references/force-fields.md 第 6.3 节
SRC_FF="${SRC_FF:-amber}"
# PBRADII: 连续介质半径集, 必须与 IGB 配套
#          1=bondi(配igb=7) 2=mbondi(igb=1) 3=mbondi2(igb=2/5) 4=mbondi3(igb=8, 默认) 7=charmm_radii(仅CHARMM)
PBRADII="${PBRADII:-4}"
# IGB: 广义 Born 模型, 1.7.0 起默认 8 (GB-Neck2); 旧版默认 5, 换默认值会导致结果不可比
IGB="${IGB:-8}"
# IDECOMP: 能量分解方式 1/2=逐残基  3/4=成对(资源消耗大)
IDECOMP="${IDECOMP:-3}"
# PRINT_RES: 官方默认并推荐 "within 6"(受体与配体间距 6A 内的残基)
#            写成 "all" 官方明确"generally not recommended"(多数残基贡献为 0)
PRINT_RES="${PRINT_RES:-within 6}"

H="${MMPBSA_HOME:-$HOME/micromamba}"
export MAMBA_ROOT_PREFIX="$H"
MM="$H/bin/micromamba"
[ -x "$MM" ] || { echo "错误: 找不到 micromamba ($MM)"; echo "先跑 scripts/setup_mmpbsa.sh"; exit 1; }

cd "$WORK" || exit 1
export GMXBIN="${GMXBIN:-/usr/local/bin}"

# ---------- 定位输入 ----------
TPR=""; XTC=""; EDR=""
for t in md.tpr prod.tpr; do [ -f "$t" ] && { TPR="$t"; break; }; done
for x in md.xtc prod.xtc traj.xtc; do [ -f "$x" ] && { XTC="$x"; break; }; done
for e in md.edr prod.edr; do [ -f "$e" ] && { EDR="$e"; break; }; done
[ -n "$TPR" ] && [ -n "$XTC" ] || { echo "错误: 找不到 tpr 或 xtc"; exit 1; }

echo "工作目录: $(pwd)"
echo "拓扑/轨迹: $TPR / $XTC"

# ---------- 1. 复合物索引: 受体组 + 配体组 + 合并组 ----------
echo
echo "[1/6] 构建索引组"
# 每次都重建: 上次残留的组名会让 dit 因重名失败, 结果组是空的(实测踩过)
rm -f idx_c.ndx idx_c2.ndx idx_r.ndx idx_l.ndx .lig_range .complex_gro

# 先确定"复合物 gro"(同时含蛋白和配体、原子数最多的那个), make_ndx 与配体定位都用它
python3 - <<'PY'
import glob
prot_res = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
            "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
best = None
for g in sorted(glob.glob("*.gro")):
    try:
        L = open(g).read().splitlines(); n = int(L[1].strip())
    except (IOError, ValueError, IndexError):
        continue
    res = [a[5:10].strip() for a in L[2:2+n]]
    lig = [i+1 for i, r in enumerate(res)
           if r not in prot_res and r not in ("SOL", "HOH", "WAT")]
    if not lig:
        continue
    if best is None or n > best[2]:
        best = (g, (min(lig), max(lig)), n)
if best is None:
    raise SystemExit("错误: 没有任何 .gro 里同时含蛋白和配体")
g, (lo, hi), n = best
print("  复合物坐标: %s (%d 原子)" % (g, n))
open(".complex_gro", "w").write(g)
open(".lig_range", "w").write("%d %d" % (lo, hi))
PY
[ -f .lig_range ] || { echo "  错误: 无法定位配体原子"; exit 1; }
COMPLEX_GRO=$(cat .complex_gro)
read LO HI < .lig_range
echo "  配体原子区间(复合物编号): $LO-$HI"

printf 'q\n' | gmx make_ndx -f "$COMPLEX_GRO" -o idx_c.ndx > /dev/null 2>&1 \
  || { echo "  错误: 无法生成默认索引"; exit 1; }

# dit 的区间是左闭右开, 要含 HI 需写 HI+1
# 坑: dit ndx_add 不会覆盖已存在的输出文件, 它会另存为 "xxx_时间戳.ndx" 并且仍返回 0,
#     结果新组静默丢失(实测踩过)。所以写之前先把目标文件删掉。
dit ndx_add -f idx_c.ndx -al "$LIG" -c "${LO}-$((HI+1))" -o idx_c2.ndx > dit1.log 2>&1 \
  || { echo "  错误: dit ndx_add 创建 $LIG 组失败"; sed 's/^/    /' dit1.log; exit 1; }
rm -f idx_c.ndx
dit ndx_add -f idx_c2.ndx -al "${PROT}_${LIG}" -c "1-$((HI+1))" -o idx_c.ndx > dit2.log 2>&1 \
  || { echo "  错误: dit ndx_add 创建合并组失败"; sed 's/^/    /' dit2.log; exit 1; }
echo "  $LIG 组: $LO-$HI ;  合并组: 1-$HI"

# ---------- 2. 受体/配体各自的局部索引(局部编号!) ----------
echo "[2/6] 生成受体/配体局部索引"
# 组名查找与 GROMACS 保持一致: 不区分大小写
# (实测: 索引里可能是 [ lig ], 传入参数却是 LIG, 大小写敏感会数到 0 个原子)
grp_size () {
  awk -v g="$1" '
    /^[[:space:]]*\[/ {
      name = tolower($0)
      sub(/^[[:space:]]*\[/, "", name)
      sub(/\][[:space:]]*$/, "", name)
      gsub(/[[:space:]]/, "", name)
      f = (name == tolower(g))
      next
    }
    f { print }
  ' idx_c.ndx | tr -s ' \t' '\n\n' | grep -c '^[0-9]\+$'
}
NPROT=$(grp_size "$PROT")
NLIG=$(grp_size "$LIG")
[ "$NPROT" -gt 0 ] || { echo "  错误: 复合物索引里找不到受体组 $PROT"; exit 1; }
[ "$NLIG" -gt 0 ]  || { echo "  错误: 复合物索引里找不到配体组 $LIG"; exit 1; }
python3 - "$PROT" "$LIG" "$NPROT" "$NLIG" <<'PY'
import sys
prot, lig, nprot, nlig = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
def wndx(fn, name, n):
    idx = list(range(1, n + 1))
    with open(fn, "w") as f:
        f.write("[ %s ]\n" % name)
        for k in range(0, len(idx), 15):
            f.write(" ".join(map(str, idx[k:k + 15])) + "\n")
wndx("idx_r.ndx", prot, nprot)
wndx("idx_l.ndx", lig, nlig)
PY
echo "  idx_r.ndx: $PROT ($NPROT 原子, 局部编号) ; idx_l.ndx: $LIG ($NLIG 原子, 局部编号)"
# ---------- 3. 拆分结构/轨迹 ----------
echo "[3/6] 拆分受体/配体结构与轨迹"
python3 - "$LO" "$HI" <<'PY'
import sys, glob, os
lo, hi = int(sys.argv[1]), int(sys.argv[2])
# 用第 1 步确定的复合物 gro(保证与配体编号一致)
if os.path.exists(".complex_gro"):
    cands = [open(".complex_gro").read().strip()]
else:
    cands = sorted(glob.glob("*.gro"))
for g in cands:
    L = open(g).read().splitlines(); n = int(L[1].strip())
    if n >= hi:
        atoms, box = L[2:2+n], L[2+n]
        open("rec.gro","w").write("receptor\n%5d\n" % (lo-1) + "\n".join(atoms[:lo-1]) + "\n" + box + "\n")
        open("lig.gro","w").write("ligand\n%5d\n" % (hi-lo+1) + "\n".join(atoms[lo-1:hi]) + "\n" + box + "\n")
        open("com.gro","w").write("complex\n%5d\n" % hi + "\n".join(atoms[:hi]) + "\n" + box + "\n")
        print("  rec.gro %d, lig.gro %d, com.gro %d 原子" % (lo-1, hi-lo+1, hi))
        break
PY
for n in com rec lig; do
  [ -f $n.gro ] && gmx editconf -f $n.gro -o $n.pdb > /dev/null 2>&1
done
printf "%s\n" "$PROT" | gmx trjconv -s "$TPR" -f "$XTC" -n idx_c.ndx -o rec.xtc > /dev/null 2>&1
printf "%s\n" "$LIG"  | gmx trjconv -s "$TPR" -f "$XTC" -n idx_c.ndx -o lig.xtc > /dev/null 2>&1
printf "%s\n" "${PROT}_${LIG}" | gmx trjconv -s "$TPR" -f "$XTC" -n idx_c.ndx -o com.xtc > /dev/null 2>&1
ls -1 rec.xtc lig.xtc com.xtc 2>/dev/null | sed 's/^/    /'

# ---------- 4. 受体/配体各自的拓扑 ----------
echo "[4/6] 生成受体/配体拓扑"
python3 - "$PROT" "$LIG" <<'PY'
import sys, glob
prot_grp, lig_grp = sys.argv[1], sys.argv[2]
for top in ("complex.top", "topol.top"):
    try: src = open(top).read().splitlines()
    except IOError: continue
    head = []
    for L in src:
        if L.strip().startswith("[ molecules ]"): break
        head.append(L)
    # 复合物拓扑: 蛋白分子 + 配体分子
    def w(fn, entries):
        open(fn,"w").write("\n".join(head + ["", "[ molecules ]", "; Compound        #mols"] + entries) + "\n")
    open("rec.top","w").write("\n".join(head + ["", "[ molecules ]", "; Compound        #mols", "Protein_chain_A     1"]) + "\n")
    open("lig.top","w").write("\n".join(head + ["", "[ molecules ]", "; Compound        #mols", "lig                 1"]) + "\n")
    print("  rec.top / lig.top 已生成")
    break
PY
EM_MDP=$(ls /usr/local/lib/python3.14/dist-packages/DuIvyTools/DuIvyTools/data/mdps/em.mdp 2>/dev/null || echo "")
[ -n "$EM_MDP" ] && {
  gmx grompp -f "$EM_MDP" -c rec.gro -p rec.top -o rec.tpr > /dev/null 2>&1
  gmx grompp -f "$EM_MDP" -c lig.gro -p lig.top -o lig.tpr > /dev/null 2>&1
}

# ---------- 5. 写 mmpbsa.in 并运行 ----------
echo "[5/6] 运行 gmx_MMPBSA"
cat > mmpbsa.in <<EOF
&general
  sys_name      = "$(basename "$WORK")",
  startframe    = $FSTART,
  endframe      = $FEND,
  interval      = 1,
  source_force_field = "$SRC_FF",
  PBRadii       = $PBRADII,
  temperature   = 298.15,
  gmx_path      = "$GMXBIN",
  solvated_trajectory = 0,
  verbose       = 2,
/

&gb
  igb         = $IGB,
  saltcon     = 0.150,
  intdiel     = 1.0,
  extdiel     = 78.5,
  surften     = 0.0072,
/

&decomp
  idecomp     = $IDECOMP,
  print_res   = "$PRINT_RES",
/
EOF
"$MM" run -n mmpbsa gmx_MMPBSA \
  -cs com.pdb -ci idx_c.ndx -cg "$PROT" "$LIG" -ct com.xtc -cp complex.top \
  -rs rec.pdb -ri idx_r.ndx -rg "$PROT" -rt rec.xtc -rp rec.top \
  -ls lig.pdb -li idx_l.ndx -lg "$LIG" -lt lig.xtc -lp lig.top \
  -i mmpbsa.in -o mmpbsa.out > mmpbsa_run.log 2>&1
RC=$?
echo "  gmx_MMPBSA exit=$RC"
if [ $RC -ne 0 ]; then
  echo "  失败原因:"
  grep -aE "MMPBSA_Error|Error:|error:" mmpbsa_run.log | head -3 | sed 's/^/    /'
  exit $RC
fi

# ---------- 6. 汇总 ----------
echo "[6/6] 结果"
python3 - <<'PY'
import re
txt = open("mmpbsa.out", errors="replace").read()

def block(tag):
    # 取小节表格的数据行(跳过标题行和分隔线)
    m = re.search(r"^" + re.escape(tag) + r".*?\n-{20,}\n(.*?)\n-{20,}", txt, re.S | re.M)
    return m.group(1) if m else ""

def total_of(b):
    # 兼容 "TOTAL" 与 "ΔTOTAL"(分解段每一行都带 Δ 前缀)
    rows = re.findall(r"^[^A-Za-z0-9]*TOTAL\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", b, re.M)
    return rows[-1] if rows else None

delta = total_of(block("Delta (Complex - Receptor - Ligand):"))
if delta is None:
    m = re.search(r"DELTA TOTAL BINDING ENERGY(.*?)\n-{10,}", txt, re.S)
    print(m.group(0).strip()[:900]) if m else print("  (未解析到结果, 详见 mmpbsa.out)")
    raise SystemExit(0)

print("  %-9s %13s %13s %13s %10s %10s" % ("", "Average", "SD", "SEM", "BlockSD", "BlockSEM"))
for name, tag in (("Complex", "Complex:"), ("Receptor", "Receptor:"), ("Ligand", "Ligand:")):
    r = total_of(block(tag))
    if r:
        print("  %-9s %13s %13s %13s %10s %10s" % (name, r[0], r[2], r[4], r[5], r[6]))
print("  %-9s %13s %13s %13s %10s %10s" % ("DELTA", delta[0], delta[2], delta[4], delta[5], delta[6]))
print()
print("  ΔG_bind = %s kcal/mol  (SEM %s, Block SD %s, Block SEM %s)" % (delta[0], delta[4], delta[5], delta[6]))
print("  注意: 'Ligand:' 段是配体自身的绝对能量, ΔG 在 'Delta (Complex - Receptor - Ligand):' 段")
print("  注意: SD(Prop.)/SEM(Prop.) 对 1e6 量级的绝对键合能会放大到无意义, 报告用 Block SD/SEM")
PY
echo
echo "完整输出: $(pwd)/mmpbsa.out"
