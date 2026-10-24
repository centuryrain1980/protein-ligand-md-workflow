#!/bin/bash
# 配体路径一条命令: SMILES -> 3D -> GAFF 参数 -> 力场合并 -> 复合物组装
# 用法: bash setup_ligand.sh <SMILES文件或mol2> <蛋白拓扑topol.top> <蛋白坐标protein.gro> <输出前缀>
set -u
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIGIN="${1:?用法: setup_ligand.sh <SMILES/mol2> <蛋白topol.top> <蛋白protein.gro> <前缀>}"
PROTTOP="${2:?}"
PROTGRO="${3:?}"
PREFIX="${4:-complex}"
FF="${PROT_FF:-amber99sb}"
WATER="${PROT_WATER:-tip3p}"
TOPDIR="${GMXDATA:-/usr/local/share/gromacs/top}"

export PATH="$HOME/.local/bin:$HOME/bin:$PATH"
cd "$(dirname "$PROTGRO")" || exit 1
WORK="$(pwd)"

echo "=== 1. SMILES -> 3D mol2 (Open Babel 3.x: -o 接格式, -O 接文件名) ==="
OB=$(command -v obabel || echo /usr/local/bin/obabel)
case "$LIGIN" in
  *.smi|*.smiles)
    "$OB" -ismi "$LIGIN" -omol2 -O"${PREFIX}_lig.mol2" --gen3d 2>&1 | tail -1
    ;;
  *.mol2)
    cp "$LIGIN" "${PREFIX}_lig.mol2"; echo "  直接使用 mol2"
    ;;
  *) echo "错误: 只支持 .smi/.smiles/.mol2"; exit 1 ;;
esac
[ -f "${PREFIX}_lig.mol2" ] || { echo "  错误: 未生成 mol2"; exit 1; }

echo "=== 2. GAFF 参数 (acpype 2026 语法) ==="
acpype -i "${PREFIX}_lig.mol2" -b lig -a gaff -o gmx -c bcc 2>&1 | tail -3
LIGDIR="$WORK/lig.acpype"
[ -f "$LIGDIR/lig_GMX.itp" ] || { echo "  错误: acpype 未生成 lig.acpype/lig_GMX.itp"; exit 1; }
echo "  产物: $(ls "$LIGDIR" | tr '\n' ' ')"

echo "=== 3. 合并 GAFF 原子类型到蛋白力场 ==="
python3 "$SKILL_DIR/scripts/merge_gaff_ff.py" \
  --lig-itp "$LIGDIR/lig_GMX.itp" \
  --ff-src "$TOPDIR/${FF}.ff" \
  --out-dir "$WORK/${FF}_gaff.ff" \
  --strip-lig-itp "$WORK/lig.itp"

echo "=== 4. 组装复合物 ==="
python3 "$SKILL_DIR/scripts/assemble_complex.py" \
  --protein "$PROTGRO" --ligand "$LIGDIR/lig_GMX.gro" \
  --topol "$PROTTOP" --lig-itp "$WORK/lig.itp" \
  --ff-dir "$WORK/${FF}_gaff.ff" --out-prefix "$PREFIX"

echo
echo "=== 完成 ==="
echo "  复合物: $WORK/${PREFIX}.gro / $WORK/${PREFIX}.top"
echo "  下一步:"
echo "    gmx editconf -f ${PREFIX}.gro -o ${PREFIX}_c.gro -c -box 4.0 4.0 4.0"
echo "    gmx solvate -cp ${PREFIX}_c.gro -cs $TOPDIR/spc216.gro -box 4.0 4.0 4.0 -o solv.gro -p ${PREFIX}.top"
