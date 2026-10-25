#!/bin/bash
# 批量执行 MD 阶段
#
# 用法:
#   bash run_simulation.sh <工作目录> [生产阶段 nsteps] [选项]
#   bash run_simulation.sh <工作目录> --plan plan.json
#
# 选项:
#   --resume        生产 MD 从检查点续跑(用 -cpi md.cpt -append, 轨迹追加不覆盖)
#   --restart       强制删除已有的生产 MD 产物(md.*)后从头跑
#   --plan <file>   用 JSON 计划驱动多阶段(见 assets/plan/*.json), 会覆盖默认的四阶段
#
# JSON 计划格式(每个 stage 一个条目):
#   name       阶段名(输出 <name>.tpr / <name>.gro / <name>.xtc)
#   mdp        mdp 模板; 相对路径先按"工作目录"找, 找不到再按"skill 目录"找
#   struct     输入坐标 (默认: 上一阶段的 .gro; 第一个阶段默认 solv_ions.gro)
#   restraint  位置限制坐标 (-r)
#   cpt        检查点 (-t)
#   nsteps     覆盖 mdp 里的 nsteps
#   continue   true 表示该阶段用 continuation = yes (接上一阶段的速度)
#   note       备注(仅供阅读)
set -u

WORK="${1:?用法: run_simulation.sh <工作目录> [生产nsteps] [--resume|--restart|--plan file]}"
if [ $# -gt 0 ]; then shift; fi
NSTEPS_MD="50000000"
if [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; then NSTEPS_MD="$1"; shift; fi

RESUME=0
RESTART=0
PLAN=""
while [ $# -gt 0 ]; do
  case "$1" in
    --resume)  RESUME=1 ;;
    --restart) RESTART=1 ;;
    --plan)    PLAN="${2:?--plan 需要一个 json 文件}"; shift ;;
    *) echo "警告: 忽略未知选项 $1"; exit 1 ;;
  esac
  shift
done

NTOMP="${NTOMP:-6}"
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MDP="$SKILL_DIR/assets/mdp"

cd "$WORK" || exit 1
echo "工作目录: $(pwd)"
echo "模板目录: $MDP"
echo "线程数  : -ntomp $NTOMP -pin on"
[ -n "$PLAN" ] && echo "阶段计划: $PLAN"
[ $RESUME -eq 1 ] && echo "续跑模式: 是 (--resume)"
[ $RESTART -eq 1 ] && echo "强制重启: 是 (--restart)"
echo

# 定位拓扑(复合物体系是 complex.top, 纯蛋白体系是 topol.top)
TOPOL=""
for t in complex.top topol.top prod.top; do
  if [ -f "$t" ]; then TOPOL="$t"; break; fi
done
[ -n "$TOPOL" ] || { echo "错误: 找不到拓扑文件 (complex.top / topol.top / prod.top)"; exit 1; }
echo "拓扑文件: $TOPOL"

INIT="${INIT:-solv_ions.gro}"
[ -f "$INIT" ] || { echo "错误: 找不到 $INIT (可设 INIT=xxx.gro)"; exit 1; }

# ---------------------------------------------------------------------------
# 阶段执行
# ---------------------------------------------------------------------------
run_stage() {
  local name="$1" mdp="$2" struct="$3" extra="$4"
  echo "---------- $name ----------"
  [ -f "$struct" ] || { echo "  错误: 找不到输入坐标 $struct"; exit 1; }
  gmx grompp -f "$mdp" -c "$struct" -p "$TOPOL" $extra -o "${name}.tpr" || {
    echo "  $name 的 grompp 失败"; exit 1; }
  local w=$(grep -ac "^WARNING" "${name}.log" 2>/dev/null || echo 0)
  gmx mdrun -s "${name}.tpr" -deffnm "$name" -ntmpi 1 -ntomp "$NTOMP" -pin on || {
    echo "  $name 的 mdrun 失败"; exit 1; }
  grep -aE "Performance:" "${name}.log" | tail -1 | sed 's/^/  /'
  echo
}

# ---------------------------------------------------------------------------
# 路径解析: mdp 优先按工作目录找, 再按 skill 目录找
# ---------------------------------------------------------------------------
resolve_mdp() {
  local p="$1"
  if [ -f "$p" ]; then echo "$p"; return 0; fi
  if [ -f "$SKILL_DIR/$p" ]; then echo "$SKILL_DIR/$p"; return 0; fi
  if [ -f "$MDP/$p" ]; then echo "$MDP/$p"; return 0; fi
  echo ""
  return 1
}

# ---------------------------------------------------------------------------
# 模式 A: JSON 计划
# ---------------------------------------------------------------------------
if [ -n "$PLAN" ]; then
  [ -f "$PLAN" ] || { echo "错误: 找不到计划文件 $PLAN"; exit 1; }
  # 让 python 把计划展平成 TSV: name|mdp|struct|restraint|cpt|nsteps|continue
  SPEC=$(python3 - "$PLAN" <<'PY'
import json, sys
try:
    with open(sys.argv[1], encoding="utf-8") as fh:
        plan = json.load(fh)
except Exception as exc:
    print("PLAN_ERROR: %s" % exc, file=sys.stderr)
    raise SystemExit(2)
if isinstance(plan, dict):
    stages = plan.get("stages", [])
else:
    stages = plan
if not stages:
    print("PLAN_ERROR: 计划里没有 stages", file=sys.stderr)
    raise SystemExit(2)
prev = None
for st in stages:
    name = st.get("name")
    if not name:
        print("PLAN_ERROR: 有 stage 缺 name", file=sys.stderr)
        raise SystemExit(2)
    struct = st.get("struct") or (f"{prev}.gro" if prev else "solv_ions.gro")
    # continuation 是 mdp 里的选项(continuation = yes/no),
    # grompp 命令行**没有** -c yes 这种用法(那个 -c 是坐标文件), 所以只能改 mdp
    cont = st.get("continue")
    if cont is None:
        cont_field = ""            # 不指定 -> 沿用 mdp 自带的设置
    elif isinstance(cont, bool):
        cont_field = "yes" if cont else "no"
    else:
        cont_field = str(cont)
    row = [name,
           st.get("mdp", ""),
           struct,
           st.get("restraint", ""),
           st.get("cpt", ""),
           str(st.get("nsteps", "")),
           cont_field]
    prev = name
    print("|".join(row))
PY
) || { echo "错误: 计划文件解析失败"; exit 1; }

  [ -n "$SPEC" ] || { echo "错误: 计划里没有有效阶段"; exit 1; }
  echo "计划阶段:"
  while IFS='|' read -r n m s r c ns cont; do
    [ -n "$n" ] || continue
    echo "  $n  mdp=$m  struct=$s  nsteps=${ns:-默认}  continuation=${cont:-跟随mdp}"
  done <<< "$SPEC"
  echo

  # 用 here-string 而不是管道: 管道里的 while 是子 shell, exit 1 拦不住后续阶段
  while IFS='|' read -r n m s r c ns cont; do
    [ -n "$n" ] || continue
    MDPP=$(resolve_mdp "$m")
    [ -n "$MDPP" ] || { echo "错误: 找不到 mdp 模板 $m"; exit 1; }
    # 阶段内覆盖 nsteps / continuation(都写进 mdp, 因为 grompp 没有对应的命令行开关)
    if [ -n "$ns" ] || [ -n "$cont" ]; then
      TMPD=".plan_$n.mdp"
      cp "$MDPP" "$TMPD"
      if [ -n "$ns" ]; then
        sed -i "s/^nsteps .*=.*/nsteps                  = $ns/" "$TMPD"
      fi
      if [ -n "$cont" ]; then
        if grep -q "^continuation" "$TMPD"; then
          sed -i "s/^continuation .*=.*/continuation            = $cont/" "$TMPD"
        else
          echo "continuation            = $cont" >> "$TMPD"
        fi
      fi
      MDPP="$TMPD"
    fi
    EXTRA=""
    [ -n "$r" ] && EXTRA="$EXTRA -r $r"
    [ -n "$c" ] && EXTRA="$EXTRA -t $c"
    run_stage "$n" "$MDPP" "$s" "$EXTRA"
  done <<< "$SPEC"
  echo "=============================================="
  echo "计划执行完成。产物在: $(pwd)"
  echo "下一步: bash <skill>/scripts/analyze.sh $(pwd)"
  exit 0
fi

# ---------------------------------------------------------------------------
# 模式 B: 默认四阶段 (EM -> NVT -> NPT -> 生产)
# ---------------------------------------------------------------------------
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