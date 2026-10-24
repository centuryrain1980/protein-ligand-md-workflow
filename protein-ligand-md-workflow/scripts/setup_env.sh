#!/bin/bash
# 环境自检与依赖安装（GROMACS 2025 + DuIvyTools 工具链）
set -u

GREEN=$'\033[32m'; YELLOW=$'\033[33m'; RED=$'\033[31m'; NC=$'\033[0m'
ok()   { printf "  ${GREEN}[OK]${NC}   %s\n" "$1"; }
warn() { printf "  ${YELLOW}[WARN]${NC} %s\n" "$1"; }
bad()  { printf "  ${RED}[FAIL]${NC} %s\n" "$1"; }

echo "=============================================="
echo " 0. 一键安装 (缺什么装什么, 可重复执行)"
echo "=============================================="
INSTALL_SH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/install_all.sh"
if [ -f "$INSTALL_SH" ]; then
  echo "  一键装齐全部环境: bash $INSTALL_SH"
  echo "  (幂等, 已装好的会跳过; 编译 GROMACS 最耗时, 可加 --without-gromacs)"
else
  echo "  未找到 install_all.sh"
fi
echo
echo "=============================================="
echo " 1. GROMACS"
echo "=============================================="
if command -v gmx >/dev/null; then
  gmx --version 2>&1 | grep -aE "GROMACS version|SIMD instructions|MPI library:" | sed 's/^/  /'
  SIMD=$(gmx --version 2>&1 | grep -a "SIMD instructions" | awk '{print $3}')
  [ "$SIMD" = "AVX2_256" ] && ok "AVX2_256 优化版 (推荐)" || warn "当前 SIMD=$SIMD"
  [ -f /usr/local/share/gromacs/top/oplsaa.ff/forcefield.itp ] \
    && ok "力场数据 /usr/local/share/gromacs/top" \
    || bad "力场数据缺失, 需重新 make install"
  [ -f /usr/local/share/gromacs/top/spc216.gro ] && ok "溶剂模板 spc216.gro" || bad "spc216.gro 缺失"
else
  bad "gmx 未安装"
fi
command -v gmx_mpi >/dev/null && ok "gmx_mpi (MPI 版)" || warn "gmx_mpi 缺失"

echo
echo "=============================================="
echo " 2. DuIvyTools (dit)"
echo "=============================================="
if command -v dit >/dev/null; then
  ok "dit 可用: $(command -v dit)"
  dit --help 2>&1 | head -1 | sed 's/^/        /'
else
  warn "dit 未安装, 安装: pip3 install --break-system-packages -i https://pypi.tuna.tsinghua.edu.cn/simple DuIvyTools"
fi

echo
echo "=============================================="
echo " 3. 化学工具"
echo "=============================================="
for t in obabel vina gnina; do
  if command -v "$t" >/dev/null; then
    ok "$t -> $(command -v $t)"
  elif [ -x "$HOME/bin/$t" ]; then
    warn "$t 在 ~/bin 但不在 PATH, 加: export PATH=\$HOME/bin:\$PATH"
  elif [ -x "/usr/local/bin/$t" ]; then
    ok "$t -> /usr/local/bin/$t"
  else
    warn "$t 缺失"
  fi
done
echo "  注: Open Babel 3.x 的语法是 -i<格式> <文件> -o<格式> -O<输出名>"

echo
echo "=============================================="
echo " 4. 配体参数化 (acpype)"
echo "=============================================="
if command -v acpype >/dev/null; then
  ok "acpype 可用: $(command -v acpype)"
elif [ -x "$HOME/.local/bin/acpype" ]; then
  warn "acpype 在 ~/.local/bin, 不在 PATH; 加: export PATH=\$HOME/.local/bin:\$PATH"
else
  warn "acpype 未安装。一键安装: bash $INSTALL_SH"
  echo "        (Ubuntu 26.04 的 rich 包由 Debian 管理, pip 需要 --break-system-packages;"
  echo "         acpype 2026 自带 antechamber, 不需要另外装 AmberTools 才能做配体参数化)"
fi

echo
echo "=============================================="
echo " 5. Python 分析库"
echo "=============================================="
python3 - <<'PY' 2>/dev/null || echo "  python3 不可用"
import importlib
for m in ["numpy","scipy","matplotlib","pandas","MDAnalysis","mdtraj"]:
    try:
        importlib.import_module(m); print(f"  \033[32m[OK]\033[0m   {m}")
    except Exception:
        print(f"  \033[33m[WARN]\033[0m {m}")
PY

echo
echo "=============================================="
echo " 6. 运行环境"
echo "=============================================="
echo "  用户: $(whoami)   (务必用普通用户, root 会导致文件属主问题)"
echo "  CPU : $(nproc) 逻辑核 / $(grep -c ^processor /proc/cpuinfo) 物理核"
echo "  内存: $(free -h | awk '/^Mem:/{print $2}')"
echo "  提示: /tmp 在 WSL 重启后会被清空, 工作目录请放在 ~ 或工作区"
echo "  提示: 混合架构 CPU, 单体系推荐 -ntomp 6 -pin on"
echo
