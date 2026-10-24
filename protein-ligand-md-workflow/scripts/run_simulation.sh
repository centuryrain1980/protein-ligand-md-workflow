#!/bin/bash
# 批量执行 EM -> NVT -> NPT -> 生产 MD
# 用法: bash run_simulation.sh <工作目录> [生产阶段 nsteps]
set -u

WORK="${1:?用法: run_simulation.sh <工作目录> [生产nsteps]}"
NSTEPS_MD="${2:-50000000}"
NTOMP="${NTOMP:-6}"
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MDP="$SKILL_DIR/assets/mdp"

cd "$WORK" || exit 1
echo "工作目录: $(pwd)"
echo "模板目录: $MDP"
echo "线程数  : -ntomp $NTOMP -pin on"
echo

# 定位拓扑与坐标
TOPOL=topol.top
[ -f "$TOPOL" ] || { echo "错误: 找不到 $TOPOL"; exit 1; }

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
sed "s/^nsteps .*=.*/nsteps                  = $NSTEPS_MD/" "$MDP/md.mdp" > md_run.mdp
run_stage md  "md_run.mdp"       "npt.gro" "-t npt.cpt"

echo "=============================================="
echo "完成。轨迹: md.xtc / md.trr / md.edr"
echo "下一步: bash <skill>/scripts/analyze.sh $(pwd)"
