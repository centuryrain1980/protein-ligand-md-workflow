#!/bin/bash
# 批量执行 EM -> NVT -> NPT -> 生产 MD
# 用法: bash run_simulation.sh <工作目录> [生产阶段 nsteps] [--resume | --restart]
#
#   --resume   生产 MD 从检查点续跑(存在 md.cpt 时用 -cpi md.cpt -append),
#              用于长模拟被中断后续跑; 其余阶段照常重跑
#   --restart  强制删除已有的生产 MD 产物(md.*)后从头跑
set -u

WORK="${1:?用法: run_simulation.sh <工作目录> [生产nsteps] [--resume|--restart]}"
NSTEPS_MD="${2:-50000000}"
shift 2 2>/dev/null || shift $#
RESUME=0
RESTART=0
for arg in "$@"; do
  case "$arg" in
    --resume)  RESUME=1 ;;
    --restart) RESTART=1 ;;
    *) echo "警告: 忽略未知选项 $arg"; exit 1 ;;
  esac
done
# 只有第 1 个参数时 shift 2 会失败, 此时没有额外选项
NTOMP="${NTOMP:-6}"
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MDP="$SKILL_DIR/assets/mdp"

cd "$WORK" || exit 1
echo "工作目录: $(pwd)"
echo "模板目录: $MDP"
echo "线程数  : -ntomp $NTOMP -pin on"
echo "生产步数: $NSTEPS_MD"
[ $RESUME -eq 1 ] && echo "续跑模式: 是 (--resume)"
[ $RESTART -eq 1 ] && echo "强制重启: 是 (--restart)"
echo

# 定位拓扑与坐标
# 复合物体系(蛋白+配体)的拓扑是 complex.top; 纯蛋白体系才是 topol.top
TOPOL=""
for t in complex.top topol.top prod.top; do
  if [ -f "$t" ]; then TOPOL="$t"; break; fi
done
[ -n "$TOPOL" ] || { echo "错误: 找不到拓扑文件 (complex.top / topol.top / prod.top)"; exit 1; }
echo "拓扑文件: $TOPOL"

run_stage() {
  local name="$1" mdp="$2" struct="$3" extra="$4"
  echo "---------- $name ----------"
  gmx grompp -f "$mdp" -c "$struct" -p "$TOPOL" $extra -o "${name}.tpr" || {
    echo "$name 的 grompp 失败"; exit 1; }
  gmx mdrun -s "${name}.tpr" -deffnm "$name" -ntmpi 1 -ntomp "$NTOMP" -pin on || {
    echo "$name 的 mdrun 失败"; exit 1; }
  grep -aE "Performance:" "${name}.log" | tail -1 | sed 's/^/  /'
  echo
}

# 初始坐标: 溶剂化+加离子后的文件
INIT="${INIT:-solv_ions.gro}"
[ -f "$INIT" ] || { echo "错误: 找不到 $INIT (可设 INIT=xxx.gro)"; exit 1; }

run_stage em  "$MDP/em.mdp"  "$INIT"  ""

# NVT/NPT 需要位置限制坐标, 用初始结构
run_stage nvt "$MDP/nvt.mdp" "em.gro" "-r $INIT"

run_stage npt "$MDP/npt.mdp" "nvt.gro" "-r $INIT -t nvt.cpt"

# 生产 MD: 无位置限制, 从 NPT 检查点续跑
if [ $RESTART -eq 1 ]; then
  echo "---------- md (强制重启, 删除旧产物) ----------"
  rm -f md.cpt md.xtc md.trr md.edr md.log md.gro md.tpr
fi

sed "s/^nsteps .*=.*/nsteps                  = $NSTEPS_MD/" "$MDP/md.mdp" > md_run.mdp

if [ $RESUME -eq 1 ] && [ -f md.cpt ]; then
  echo "---------- md (从检查点续跑) ----------"
  if [ -f md.tpr ]; then
    old=$(grep -aE "^ *nsteps " md.tpr | head -1 | tr -s ' ' | cut -d' ' -f3)
    if [ -n "$old" ] && [ "$old" != "$NSTEPS_MD" ]; then
      echo "  提示: 已有 md.tpr 的 nsteps=$old, 本次设为 $NSTEPS_MD"
      echo "        (续跑时 GROMACS 会以 md.tpr 里的设置为准, 想改总步数请先 --restart)"
      NSTEPS_MD="$old"
      sed "s/^nsteps .*=.*/nsteps                  = $NSTEPS_MD/" "$MDP/md.mdp" > md_run.mdp
    fi
  fi
  gmx grompp -f md_run.mdp -c npt.gro -t npt.cpt -p "$TOPOL" -o md.tpr || {
    echo "md 的 grompp 失败"; exit 1; }
  # -cpi 从检查点续跑, -append 把新帧追加到原轨迹(否则会覆盖)
  gmx mdrun -s md.tpr -deffnm md -cpi md.cpt -append -ntmpi 1 -ntomp "$NTOMP" -pin on || {
    echo "md 的 mdrun 失败"; exit 1; }
  grep -aE "Performance:" md.log | tail -1 | sed 's/^/  /'
  echo
else
  run_stage md  "md_run.mdp"       "npt.gro" "-t npt.cpt"
fi

echo "=============================================="
echo "完成。轨迹: md.xtc / md.trr / md.edr"
echo "下一步: bash <skill>/scripts/analyze.sh $(pwd)"
echo "或用统一入口: mdflow.py analyze --workdir $(pwd)"