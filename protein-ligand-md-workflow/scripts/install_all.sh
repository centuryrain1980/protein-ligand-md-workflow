#!/bin/bash
# ==============================================================================
# 一键安装 protein-ligand-md-workflow 的全部运行环境
#
#   bash install_all.sh [选项]
#
#   --with-gromacs      编译安装 GROMACS 2025(默认开; 最耗时的一步, 约 20-60 分钟)
#   --without-gromacs   跳过 GROMACS(若系统已装好 gmx 则自动跳过)
#   --with-docking      安装 Vina + gnina(默认开, 需下载约 1.5 GB)
#   --without-docking   跳过对接工具
#   --without-mmpbsa    跳过 MM/PBSA 环境(micromamba + AmberTools, 约 2 GB)
#   --with-charmm36     额外下载 MacKerell 的 CHARMM36 力场文件(14 MB)
#                        注意: 该 ff 与 GROMACS 2025 的 pdb2gmx 不兼容,
#                        GROMACS 2025 上做 CHARMM 请用自带的 charmm27.ff
#   --mirror <url>      pip 镜像(默认清华源)
#   -h | --help         显示本帮助
#
# 幂等: 已装好的组件会自动跳过, 可以重复执行。
# 依赖: 需要 sudo 权限(装系统包和 GROMACS 到 /usr/local)。
# ==============================================================================
set -u

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

GROMACS_VERSION="${GROMACS_VERSION:-2025.4}"
PIP_INDEX_URL="${PIP_INDEX_URL:-https://pypi.tuna.tsinghua.edu.cn/simple}"
MAMBA_TARBALL_URL="${MAMBA_TARBALL_URL:-https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/micromamba-2.8.1-1.tar.bz2}"
VINA_URL="${VINA_URL:-https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64}"
GNINA_URL="${GNINA_URL:-https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2}"

WITH_GROMACS=1
WITH_DOCKING=1
WITH_MMPBSA=1
WITH_CHARMM36=0
BUILD_JOBS="${BUILD_JOBS:-$(nproc 2>/dev/null || echo 4)}"

# ---------- 输出helper ----------
if [ -t 1 ]; then G=$'\033[32m'; Y=$'\033[33m'; R=$'\033[31m'; B=$'\033[1m'; N=$'\033[0m'
else G=""; Y=""; R=""; B=""; N=""; fi
log()  { printf "  ${G}[OK]${N}   %s\n" "$1"; }
warn() { printf "  ${Y}[WARN]${N} %s\n" "$1"; }
err()  { printf "  ${R}[FAIL]${N} %s\n" "$1"; }
step() { printf "\n${B}==> %s${N}\n" "$1"; }
die()  { err "$1"; exit 1; }

usage() { sed -n '2,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

while [ $# -gt 0 ]; do
  case "$1" in
    --with-gromacs)    WITH_GROMACS=1 ;;
    --without-gromacs) WITH_GROMACS=0 ;;
    --with-docking)    WITH_DOCKING=1 ;;
    --without-docking) WITH_DOCKING=0 ;;
    --without-mmpbsa)  WITH_MMPBSA=0 ;;
    --with-charmm36)   WITH_CHARMM36=1 ;;
    --mirror)          PIP_INDEX_URL="${2:?--mirror 需要一个 URL}"; shift ;;
    -h|--help)         usage ;;
    *) die "未知选项: $1 (用 --help 看用法)" ;;
  esac
  shift
done

echo "${B}=========================================================${N}"
echo "${B} protein-ligand-md-workflow 环境安装${N}"
echo "${B}=========================================================${N}"
echo "  GROMACS 编译 : $([ $WITH_GROMACS -eq 1 ] && echo 开 || echo 关)"
echo "  对接工具     : $([ $WITH_DOCKING -eq 1 ] && echo 开 || echo 关)"
echo "  MM/PBSA 环境 : $([ $WITH_MMPBSA -eq 1 ] && echo 开 || echo 关)"
echo "  CHARMM36 ff   : $([ $WITH_CHARMM36 -eq 1 ] && echo 开 || echo 关)"
echo "  pip 镜像     : $PIP_INDEX_URL"
echo "  并行编译数   : $BUILD_JOBS"

# ---------- sudo ----------
SUDO=""
if [ "$(id -u)" -eq 0 ]; then
  SUDO=""
elif command -v sudo >/dev/null 2>&1; then
  # 先认证一次把凭证缓存住; 之后一律用 sudo -n, 这样在无终端(自动化/CI)环境下
  # 会立刻失败而不是卡在"请输入密码"上
  # sudo 默认从 /dev/tty 读密码, 重定向 stdin 拦不住, 所以按有没有终端分别处理:
  #   有终端 -> 交互输一次密码并缓存凭证
  #   无终端 -> 一律用 sudo -n, 需要密码就立刻失败, 绝不挂起
  if [ -t 0 ] && [ -t 1 ]; then
    if ! timeout 120 sudo -v; then
      die "sudo 认证失败, 无法继续安装系统包。"
    fi
  elif ! sudo -n -v 2>/dev/null; then
    warn "当前不是交互式终端且 sudo 需要密码, 特权步骤会失败并跳过。"
    warn "  · 在交互式终端里重跑本脚本即可正常完成;"
    warn "  · 或让管理员先装好系统包与 GROMACS, 再用 --without-gromacs 运行。"
  fi
  SUDO="sudo -n"
else
  die "需要 sudo 权限来装系统包和 GROMACS, 但找不到 sudo 命令"
fi

run_root() { if [ -n "$SUDO" ]; then $SUDO "$@" </dev/null; else "$@" </dev/null; fi; }

# ---------- 0. 前置检查 ----------
step "0. 前置检查"
if ! command -v apt-get >/dev/null 2>&1; then
  die "只支持 Debian/Ubuntu 系(检测不到 apt-get)。
     Rocky/CentOS 请把 setup_env.sh 里的 apt 换成 dnf, 其余步骤一致。"
fi
log "包管理器: apt-get"

if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
  warn "没有 curl/wget, 尝试用 python3 下载"
fi
if ! curl -s -o /dev/null --max-time 20 https://pypi.org/simple/ 2>/dev/null; then
  warn "连不上 pypi.org (可能需要代理)。如有代理请先 export HTTPS_PROXY=..."
fi

# ---------- 1. 系统依赖 ----------
step "1. 系统依赖 (apt)"
if [ "${SKIP_APT:-0}" = "1" ]; then
  warn "SKIP_APT=1, 跳过 apt 步骤"
else
SYS_PKGS="build-essential cmake make wget curl bzip2 ca-certificates python3-pip python3-venv libfftw3-dev libomp-dev"
MISSING=""
for p in $SYS_PKGS; do
  dpkg -s "$p" >/dev/null 2>&1 || MISSING="$MISSING $p"
done
if [ -z "$MISSING" ]; then
  log "系统依赖已齐 (跳过 apt 安装)"
else
  echo "  缺少:$MISSING"
  run_root apt-get update -qq || die "apt-get update 失败"
  # shellcheck disable=SC2086
  run_root apt-get install -y -qq $MISSING || die "apt-get install 失败"
  log "系统依赖安装完成"
fi
if dpkg -s openbabel >/dev/null 2>&1 || command -v obabel >/dev/null 2>&1; then
  log "Open Babel 就绪 ($(command -v obabel || echo 已装))"
else
  warn "Open Babel 未装, 对接时需要 (apt-get install openbabel)"
fi
fi  # SKIP_APT

# ---------- 2. GROMACS ----------
step "2. GROMACS $GROMACS_VERSION"
gmx_ok() {
  command -v gmx >/dev/null 2>&1 || return 1
  gmx --version 2>/dev/null | grep -q "GROMACS version: *$GROMACS_VERSION" || return 1
  [ -f /usr/local/share/gromacs/top/spc216.gro ] || [ -f "$(gmx --version 2>/dev/null | grep 'Data prefix' | awk '{print $3}')/share/gromacs/top/spc216.gro" ] || return 1
  return 0
}
if gmx_ok; then
  log "GROMACS $GROMACS_VERSION 已就绪 ($(command -v gmx))"
  gmx --version 2>/dev/null | grep -E "SIMD instructions|MPI library" | sed 's/^/        /'
elif [ "$WITH_GROMACS" -eq 0 ]; then
  warn "按要求跳过 GROMACS 安装 (当前: $(command -v gmx || echo 未安装))"
else
  SIMD="AVX2_256"
  if ! grep -qw avx2 /proc/cpuinfo 2>/dev/null; then
    SIMD="AVX2"
    grep -qw avx /proc/cpuinfo 2>/dev/null || SIMD="SSE4.2"
  fi
  echo "  目标 SIMD: $SIMD"
  BUILD_ROOT="$HOME/build"
  SRC="$BUILD_ROOT/gromacs-$GROMACS_VERSION"
  mkdir -p "$BUILD_ROOT"
  if [ ! -d "$SRC" ]; then
    echo "  下载源码包 (~90 MB)"
    DL="$BUILD_ROOT/gromacs-$GROMACS_VERSION.tar.gz"
    if [ ! -f "$DL" ]; then
      curl -fL --retry 3 --connect-timeout 30 -o "$DL" \
        "https://ftp.gromacs.org/pub/gromacs/gromacs-$GROMACS_VERSION.tar.gz" \
        || die "下载 GROMACS 源码失败"
    fi
    tar xzf "$DL" -C "$BUILD_ROOT" || die "解压失败"
  fi
  [ -d "$SRC" ] || die "源码目录不存在: $SRC"
  mkdir -p "$SRC/build" && cd "$SRC/build" || die "无法进入 build 目录"
  echo "  cmake 配置 (含 GPU=OFF, 自带 FFTW) ..."
  cmake .. \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DGMX_SIMD="$SIMD" \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_GPU=OFF \
    -DGMXAPI=OFF \
    -DREGRESSIONTEST_DOWNLOAD=OFF \
    -DGMX_BUILD_UNITTESTS=OFF \
    > cmake.log 2>&1 || { tail -25 cmake.log; die "cmake 配置失败"; }
  echo "  编译中 (make -j$BUILD_JOBS, 约 20-60 分钟) ..."
  make -j"$BUILD_JOBS" > make.log 2>&1 || { tail -25 make.log; die "编译失败"; }
  run_root make install > install.log 2>&1 || {
    tail -25 install.log
    die "make install 失败(通常是没有 sudo 权限装到 /usr/local)。
       可改为装到用户目录: cmake .. -DCMAKE_INSTALL_PREFIX=\$HOME/.local, 然后把 gmx 加进 PATH。"
  }
  hash -
  log "GROMACS 安装完成: $(command -v gmx || echo /usr/local/bin/gmx)"
  gmx --version 2>/dev/null | grep -E "GROMACS version|SIMD instructions" | sed 's/^/        /'
fi

# ---------- 3. Python 包 ----------
step "3. Python 分析与参数化工具"
PIP_PKGS="DuIvyTools acpype MDAnalysis mdtraj meeko rdkit gemmi ProDy numpy scipy pandas matplotlib"
# rich 在 Ubuntu 26.04 由 Debian 托管, 直接装会被拒; --ignore-installed 强制覆盖
pip_install() {
  python3 -m pip install --break-system-packages --ignore-installed rich \
    -i "$PIP_INDEX_URL" "$@" 2>&1 | tail -3
}
if python3 -c "import acpype, MDAnalysis, mdtraj, meeko" >/dev/null 2>&1; then
  log "Python 包已就绪 (acpype / MDAnalysis / mdtraj / meeko 均可导入)"
else
  echo "  安装: $PIP_PKGS"
  pip_install $PIP_PKGS || die "pip 安装失败 (可换镜像: --mirror https://pypi.org/simple)"
  log "Python 包安装完成"
fi
python3 - <<'PY' 2>/dev/null || true
import importlib
missing = []
for m in ["acpype", "MDAnalysis", "mdtraj", "meeko", "rdkit", "gemmi", "prody",
          "numpy", "scipy", "pandas", "matplotlib"]:
    try:
        importlib.import_module(m)
    except Exception:
        missing.append(m)
print("  仍缺失:", " ".join(missing) if missing else "无")
PY
command -v dit >/dev/null 2>&1 && log "dit 可用: $(command -v dit)" || warn "dit 未进入 PATH"

# ---------- 4. MM/PBSA 环境 ----------
step "4. MM/PBSA 环境 (micromamba + AmberTools + gmx_MMPBSA)"
if [ "$WITH_MMPBSA" -eq 0 ]; then
  warn "按要求跳过 MM/PBSA 环境"
elif [ -x "${MMPBSA_HOME:-$HOME/micromamba}/bin/micromamba" ]; then
  log "micromamba 已存在, 继续检查环境"
  bash "$SKILL_DIR/scripts/setup_mmpbsa.sh" || warn "MM/PBSA 环境检查有告警"
else
  bash "$SKILL_DIR/scripts/setup_mmpbsa.sh" || warn "MM/PBSA 环境安装有告警(可稍后重跑本脚本)"
fi

# ---------- 5. 对接工具 ----------
step "5. 对接工具 (Vina + gnina)"
if [ "$WITH_DOCKING" -eq 0 ]; then
  warn "按要求跳过对接工具"
else
  mkdir -p "$HOME/apps/vina" "$HOME/apps/gnina" "$HOME/bin"
  # --- Vina ---
  if [ -x "$HOME/apps/vina/vina-bin" ]; then
    log "Vina 已就绪"
  else
    echo "  下载 Vina 1.2.7 ..."
    if curl -fL --retry 3 --connect-timeout 30 -o "$HOME/apps/vina/vina-bin" "$VINA_URL"; then
      chmod +x "$HOME/apps/vina/vina-bin"
      ln -sf "$HOME/apps/vina/vina-bin" "$HOME/bin/vina"
      log "Vina 安装完成"
    else
      warn "Vina 下载失败, 可手动下载: $VINA_URL"
    fi
  fi
  # --- gnina (CPU 版, 不依赖 CUDA) ---
  GNINA_BIN=""
  for c in "$HOME/apps/gnina/gnina_v1.3.3" "$HOME/apps/gnina/gnina_v1.3.2" "$HOME/apps/gnina/gnina_v1.1"; do
    [ -x "$c" ] && { GNINA_BIN="$c"; break; }
  done
  if [ -n "$GNINA_BIN" ]; then
    log "gnina 已就绪: $GNINA_BIN"
    if [ ! -e "$HOME/bin/gnina" ] && [ ! -L "$HOME/bin/gnina" ]; then
      ln -sf "$GNINA_BIN" "$HOME/bin/gnina"
    fi
  else
    echo "  下载 gnina 1.3.2 (CPU 版, 约 1.4 GB, 请耐心) ..."
    if curl -fL --retry 3 --connect-timeout 30 -o "$HOME/apps/gnina/gnina_v1.3.2" "$GNINA_URL"; then
      chmod +x "$HOME/apps/gnina/gnina_v1.3.2"
      ln -sf "$HOME/apps/gnina/gnina_v1.3.2" "$HOME/bin/gnina"
      log "gnina 安装完成"
    else
      warn "gnina 下载失败, 可手动下载: $GNINA_URL"
    fi
  fi
fi

# ---------- 5b. CHARMM36 力场文件(可选) ----------
if [ "$WITH_CHARMM36" -eq 1 ]; then
  step "5b. CHARMM36 力场文件 (MacKerell)"
  CH36_URL="${CH36_URL:-https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_feb26.tgz}"
  CDIR="$HOME/apps/charmm36"
  mkdir -p "$CDIR" && cd "$CDIR" || die "无法创建 $CDIR"
  if [ -d "$CDIR/charmm36.ff" ]; then
    log "CHARMM36 力场文件已存在"
  else
    echo "  下载 toppar (~14 MB) ..."
    if curl -fL --retry 3 --connect-timeout 30 -o toppar_c36.tgz "$CH36_URL"; then
      tar xzf toppar_c36.tgz 2>/dev/null || true
      FFGZ=$(find . -path "*gromacs*" -name "charmm36-*.ff.tgz" ! -name "*ljpme*" | head -1)
      if [ -z "$FFGZ" ]; then
        warn "包里没找到 charmm36-*.ff.tgz, 请到 $CH36_URL 手工下载"
      else
        rm -rf _ff && mkdir _ff
        tar xzf "$FFGZ" -C _ff
        SRCDIR=$(find _ff -maxdepth 1 -type d -name "*.ff" | head -1)
        if [ -z "$SRCDIR" ]; then
          warn "解压后没找到 *.ff 目录"
        else
          rm -rf "$CDIR/charmm36.ff"
          cp -r "$SRCDIR" "$CDIR/charmm36.ff"
          log "CHARMM36 力场文件已装到 $CDIR/charmm36.ff"
          echo "        免 sudo 使用: export GMXLIB=/usr/local/share/gromacs/top:$CDIR"
          echo "        或装到系统库: sudo cp -r $CDIR/charmm36.ff /usr/local/share/gromacs/top/"
        fi
      fi
    else
      warn "CHARMM36 下载失败: $CH36_URL"
    fi
  fi
  warn "注意: MacKerell 的 charmm36.ff 与 GROMACS 2025 的 pdb2gmx 不兼容"
  warn "      (实测报 'atom C1 not found in buiding block 1MET while combining tdb and rtp')。"
  warn "      GROMACS 2025 上做 CHARMM 请用自带的 charmm27.ff + -water tips3p,"
  warn "      或用 CHARMM-GUI 建模, 或把 GROMACS 降到 2023 及更早。"
  cd - >/dev/null || true
fi

# ---------- 6. PATH 持久化 ----------
step "6. PATH 持久化 (~/.bashrc)"
BASHRC="$HOME/.bashrc"
touch "$BASHRC"
for line in 'export PATH="$HOME/bin:$HOME/.local/bin:$PATH"' 'export PATH="$HOME/micromamba/bin:$PATH"'; do
  if ! grep -qF "$line" "$BASHRC"; then
    printf '\n# protein-ligand-md-workflow\n%s\n' "$line" >> "$BASHRC"
    log "已加入: $line"
  else
    log "已存在: $line"
  fi
done
log "重开终端或执行: source ~/.bashrc"

# ---------- 7. 自检 ----------
step "7. 环境自检"
bash "$SKILL_DIR/scripts/setup_env.sh"

echo
echo "${B}=========================================================${N}"
echo "${B} 安装流程结束${N}"
echo "${B}=========================================================${N}"
echo "  下一步:"
echo "    1) 重开终端 (或 source ~/.bashrc) 让 PATH 生效"
echo "    2) 跑一个测试体系:"
echo "         cd ~/projects/demo && mkdir -p mdp && cd mdp"
echo "         gmx pdb2gmx -f receptor.pdb -o protein.gro -i posre.itp -p topol.top \\"
echo "                     -ff amber99sb -water tip3p -ignh < /dev/null"
echo "    3) 详细说明见 README.md"