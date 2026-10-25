#!/bin/bash
# 轨迹分析全套: 预处理 -> 常规量 -> 二级结构 -> Ramachandran -> 能量分解 -> 出图
# 用法: bash analyze.sh <工作目录> [蛋白组名] [配体组名]
set -u
export MPLBACKEND=Agg

WORK="${1:?用法: analyze.sh <工作目录> [蛋白组] [配体组]}"
PROT="${2:-Protein}"
LIG="${3:-LIG}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORK" || exit 1
mkdir -p analysis
echo "工作目录: $(pwd)   蛋白组: $PROT   配体组: $LIG"

TPR=""; XTC=""; EDR=""
for t in md.tpr topol.tpr prod.tpr; do [ -f "$t" ] && { TPR="$t"; break; }; done
for x in md.xtc prod.xtc traj.xtc; do [ -f "$x" ] && { XTC="$x"; break; }; done
for e in md.edr prod.edr; do [ -f "$e" ] && { EDR="$e"; break; }; done

if [ -z "$TPR" ] || [ -z "$XTC" ]; then
  echo "错误: 找不到 tpr 或 xtc (期望 md.tpr/md.xtc 或 prod.tpr/prod.xtc)"
  exit 1
fi
echo "拓扑: $TPR   轨迹: $XTC   能量: ${EDR:-无}"
echo

step() { echo "---------- $1 ----------"; }

# 1) 轨迹预处理 —— GROMACS 2025 要求 -pbc 与 -fit 分两步, 且 trjconv 需要多次组选择
step "轨迹预处理 (先 -pbc whole, 再 -center -fit)"
printf "%s\n%s\n" "$PROT" "$PROT" | gmx trjconv -s "$TPR" -f "$XTC" \
    -o analysis/whole.xtc -pbc whole > analysis/trjconv1.log 2>&1
if [ -f analysis/whole.xtc ]; then
  printf "%s\n%s\n%s\n" "$PROT" "$PROT" "$PROT" | gmx trjconv -s "$TPR" -f analysis/whole.xtc \
      -o analysis/traj_fit.xtc -center -fit rot+trans > analysis/trjconv2.log 2>&1
fi
if [ -f analysis/traj_fit.xtc ]; then
  echo "  完成 (去 PBC 伪影 + 按 $PROT 拟合)"
  XT=analysis/traj_fit.xtc
else
  echo "  预处理未成功, 退回原始轨迹"
  echo "  诊断: $(grep -aA3 'Fatal error' analysis/trjconv1.log analysis/trjconv2.log 2>/dev/null | head -4 | tr '\n' ' ')"
  XT="$XTC"
fi

# 2) 常规分析
step "RMSD"
printf "%s\n%s\n" "$PROT" "$PROT" | gmx rms -s "$TPR" -f "$XT" -o analysis/rmsd.xvg -tu ps > analysis/rms.log 2>&1
echo "  行数: $(grep -vc '^@' analysis/rmsd.xvg 2>/dev/null)"

step "RMSF (按残基)"
printf "%s\n" "$PROT" | gmx rmsf -s "$TPR" -f "$XT" -o analysis/rmsf.xvg -res > analysis/rmsf.log 2>&1
echo "  行数: $(grep -vc '^@' analysis/rmsf.xvg 2>/dev/null)"

step "回转半径"
printf "%s\n" "$PROT" | gmx gyrate -s "$TPR" -f "$XT" -o analysis/gyrate.xvg > analysis/gyrate.log 2>&1
echo "  行数: $(grep -vc '^@' analysis/gyrate.xvg 2>/dev/null)"

step "溶剂可及面积"
printf "%s\n" "$PROT" | gmx sasa -s "$TPR" -f "$XT" -o analysis/sasa.xvg > analysis/sasa.log 2>&1
echo "  行数: $(grep -vc '^@' analysis/sasa.xvg 2>/dev/null)"

step "氢键"
# GROMACS 2025: 新版 gmx hbond 出统计, 旧版 hbond-legacy 出氢键图谱 xpm
printf "%s\n%s\n" "$PROT" "$PROT" | gmx hbond -s "$TPR" -f "$XT" \
    -num analysis/hbond_num.xvg -dist analysis/hbond_dist.xvg > analysis/hbond.log 2>&1
[ -f analysis/hbond_num.xvg ] && echo "  统计完成" || echo "  新版 hbond 失败"
printf "%s\n%s\n" "$PROT" "$PROT" | gmx hbond-legacy -s "$TPR" -f "$XT" \
    -num analysis/hbond_num2.xvg -hbm analysis/hbond.xpm > analysis/hbond_legacy.log 2>&1
[ -f analysis/hbond.xpm ] && echo "  图谱完成 (hbond.xpm)" || echo "  hbond-legacy 不可用"

# 3) 二级结构
step "二级结构 (gmx dssp -> dit dssp)"
printf "%s\n" "$PROT" | gmx dssp -s "$TPR" -f "$XT" -o analysis/dssp.dat > analysis/dssp.log 2>&1
if [ -f analysis/dssp.dat ] && [ -s analysis/dssp.dat ]; then
  dit dssp -f analysis/dssp.dat -o analysis/dssp_out > analysis/dit_dssp.log 2>&1
  ls -1 analysis/dssp_out*.xpm analysis/dssp_out*.xvg 2>/dev/null | sed 's/^/  /'
else
  echo "  gmx dssp 失败或无输出"
fi

# 4) Ramachandran
step "Ramachandran"
printf "%s\n%s\n" "$PROT" "$PROT" | gmx rama -s "$TPR" -f "$XT" -o analysis/rama.xvg > analysis/rama.log 2>&1
echo "  行数: $(grep -vc '^@' analysis/rama.xvg 2>/dev/null)"

# 5) 能量分解
if [ -n "$EDR" ]; then
  step "能量分解 (gmx energy 项 7,8,9,10)"
  printf '7\n8\n9\n10\n' | gmx energy -f "$EDR" -o analysis/energy_sys.xvg > analysis/energy.log 2>&1
  echo "  列名: $(grep -a '^@ s0' analysis/energy_sys.xvg 2>/dev/null | head -1)"

  step "组间相互作用能 (需 mdp 里设过 energygrps)"
  printf 'q\n' | gmx energy -f "$EDR" -o /dev/null 2>&1 | grep -aE "Coul-SR:|LJ-SR:" | head -8 | sed 's/^/  /'

  if [ -f analysis/energy_sys.xvg ]; then
    cp analysis/energy_sys.xvg analysis/prolig.xvg
    cp analysis/energy_sys.xvg analysis/pro.xvg
    cp analysis/energy_sys.xvg analysis/lig.xvg
    dit xvg_energy_compute -f analysis/prolig.xvg analysis/pro.xvg analysis/lig.xvg \
        -o analysis/binding.xvg > analysis/dit_energy.log 2>&1
    [ -f analysis/binding.xvg ] && echo "  binding.xvg 已生成" || echo "  dit xvg_energy_compute 失败"
  fi
else
  echo "未找到 .edr, 跳过能量分析"
fi

# 6) 出图
step "绘图"
plot() {
  local src="$1" out="$2"
  [ -f "$src" ] || { echo "  跳过 (缺 $src)"; return; }
  dit xvg_show -f "$src" -o "$out" > analysis/plot.log 2>&1
  [ -f "$out" ] && echo "  $out" || echo "  $out 生成失败"
}
plot analysis/rmsd.xvg      analysis/rmsd.png
plot analysis/rmsf.xvg      analysis/rmsf.png
plot analysis/gyrate.xvg    analysis/gyrate.png
plot analysis/sasa.xvg      analysis/sasa.png
# Ramachandran 必须画散点: gmx rama 输出的是逐点 (phi,psi),
# dit xvg_show 只会按顺序连成折线, 出来是一团乱麻 -> 用专门的 plot_rama.py
if [ -f analysis/rama.xvg ]; then
  if [ -f "$SCRIPT_DIR/plot_rama.py" ]; then
    python3 "$SCRIPT_DIR/plot_rama.py" analysis/rama.xvg \
        -o analysis/rama.png --stat analysis/rama_stat.txt 2>&1 | sed "s/^/  /"
  else
    echo "  警告: 找不到 plot_rama.py, Ramachandran 不会出图"
  fi
fi
plot analysis/binding.xvg   analysis/binding.png
plot analysis/hbond_num.xvg analysis/hbond.png
# 体系能量分量 (范德华/离散校正/库仑/长程静电): 稳定性检查用
plot analysis/energy_sys.xvg analysis/energy.png

if [ -f analysis/dssp_out.xpm ]; then
  dit xpm_show -f analysis/dssp_out.xpm -o analysis/dssp.png > analysis/plot.log 2>&1
  [ -f analysis/dssp.png ] && echo "  analysis/dssp.png"
fi

echo
echo "=============================================="
echo "分析完成, 结果在: $(pwd)/analysis"
ls -1 analysis/*.png 2>/dev/null | sed 's/^/  /'
