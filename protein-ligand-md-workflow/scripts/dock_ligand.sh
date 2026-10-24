#!/bin/bash
# 分子对接 (gnina, 带 CNN 重打分)
# 用法: bash dock_ligand.sh <受体PDB> <配体mol2/sdf> [盒子中心X Y Z] [盒子边长]
# 说明:
#   - gnina 直接吃 PDB 受体(内部自动转 PDBQT), 不需要 meeko
#   - 配体需先转成 PDBQT (obabel)
#   - 不给盒子中心时用受体几何中心(盲对接), 建议边长 30
set -u
export PATH="$HOME/bin:$HOME/.local/bin:$PATH"

RECP="${1:?用法: dock_ligand.sh <受体PDB> <配体mol2/sdf> [cx cy cz] [size]}"
LIGP="${2:?}"
CENTER=("${@:3}")
SIZE="${CENTER[3]:-22}"
CENTER=("${CENTER[@]:0:3}")

OUT="dock_$(basename "${LIGP%.*}")"
OB=$(command -v obabel || echo /usr/local/bin/obabel)
[ -f "$RECP" ] || { echo "错误: 受体不存在 $RECP"; exit 1; }

echo "=== 受体: $RECP ($(grep -c '^ATOM' "$RECP") 原子) ==="
echo "=== 配体: $LIGP ==="

# ---- 1. 配体 -> PDBQT ----
case "$LIGP" in
  *.pdbqt) cp "$LIGP" "${OUT}_lig.pdbqt";;
  *)
    # Open Babel 3.x: -i<格式> <输入文件> -o<格式> -O<输出文件>
    case "$LIGP" in
      *.mol2)  FMT="-imol2" ;;
      *.sdf)   FMT="-isdf" ;;
      *.pdb)   FMT="-ipdb" ;;
      *.smi|*.smiles) FMT="-ismi" ;;
      *) FMT="-imol2" ;;
    esac
    "$OB" "$FMT" "$LIGP" -opdbqt -O"${OUT}_lig.pdbqt" 2>&1 | tail -1
    ;;
esac
[ -s "${OUT}_lig.pdbqt" ] || { echo "错误: 配体 PDBQT 转换失败"; exit 1; }
echo "  配体 PDBQT: $(grep -c '^ATOM' "${OUT}_lig.pdbqt") 原子"

# ---- 2. 盒子中心 ----
if [ ${#CENTER[@]} -lt 3 ] || [ -z "${CENTER[0]}" ]; then
  read CX CY CZ < <(python3 - "$RECP" <<'PY'
import sys
c=[(float(l[30:38]),float(l[38:46]),float(l[46:54])) for l in open(sys.argv[1]) if l.startswith("ATOM")]
n=len(c)
print(sum(x for x,_,_ in c)/n, sum(y for _,y,_ in c)/n, sum(z for _,_,z in c)/n)
PY
)
  SIZE=30
  echo "  未指定盒子: 用受体几何中心 (盲对接), 边长 $SIZE"
else
  CX="${CENTER[0]}"; CY="${CENTER[1]}"; CZ="${CENTER[2]}"
  echo "  盒子中心: ($CX, $CY, $CZ)  边长: $SIZE"
fi

# ---- 3. gnina 对接 ----
echo
echo "=== gnina 对接 (CPU) ==="
EXH="${EXHAUSTIVENESS:-8}"
NMODES="${NUM_MODES:-5}"
gnina --receptor "$RECP" --ligand "${OUT}_lig.pdbqt" \
      --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
      --size_x "$SIZE" --size_y "$SIZE" --size_z "$SIZE" \
      --exhaustiveness "$EXH" --num_modes "$NMODES" \
      --no_gpu --out "${OUT}_docked.sdf" --log "${OUT}_gnina.log" 2>&1 | tail -12
RC=$?
echo "  exit=$RC"
[ -f "${OUT}_docked.sdf" ] || { echo "错误: 未生成对接结果"; exit 1; }

# ---- 4. 汇总 ----
echo
echo "=== 结果 (按 CNN pose score 排序) ==="
python3 - "${OUT}_gnina.log" <<'PY'
import sys, re
txt = open(sys.argv[1], errors="replace").read()
rows = re.findall(r"^\s*(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$", txt, re.M)
if rows:
    print("  mode  affinity   intramol  CNN_pose  CNN_aff")
    for r in rows:
        print("  %4s  %9s  %9s  %9s  %9s" % tuple(r))
    best = max(rows, key=lambda r: float(r[4]))
    print("\n  CNN 最优: mode %s  affinity %s kcal/mol  CNN pose %s" % (best[0], best[1], best[3]))
else:
    print("  未解析到结果表, 见 %s" % sys.argv[1])
PY
echo
echo "产物: ${OUT}_docked.sdf  (多构象, 用 \$\$\$\$ 分隔)"
echo "下一步: 取最优构象 -> acpype 生成 GAFF 参数 -> scripts/setup_ligand.sh 组装"
