#!/bin/bash
# 安装 gmx_MMPBSA 运行环境(Python 3.12 + AmberTools + gmx_MMPBSA)
# 用法: bash setup_mmpbsa.sh
set -u
H="${MMPBSA_HOME:-$HOME/micromamba}"
export MAMBA_ROOT_PREFIX="$H"
MM="$H/bin/micromamba"

echo "=== 1. micromamba ==="
MM_URL="${MM_URL:-https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/micromamba-2.8.1-1.tar.bz2}"
if [ ! -x "$MM" ]; then
  mkdir -p "$H"
  echo "  准备 micromamba ..."
  TAR=$(ls -1 "$HOME"/../*/work/ref/micromamba.tar.bz2 2>/dev/null | head -1)
  [ -z "$TAR" ] && TAR=$(find /mnt -maxdepth 6 -name "micromamba.tar.bz2" 2>/dev/null | head -1)
  if [ -n "$TAR" ]; then
    echo "  用本地包: $TAR"
  else
    echo "  本地没有, 从镜像下载: $MM_URL"
    TAR="$H/micromamba.tar.bz2"
    curl -fL --retry 3 --connect-timeout 30 -o "$TAR" "$MM_URL" \
      || { echo "  错误: 下载 micromamba 失败, 可手动下载后放到任意目录再重跑"; exit 1; }
  fi
  echo "  解压 ..."
  tar xjf "$TAR" -C "$H" bin/micromamba || { echo "  错误: 解压失败"; exit 1; }
  chmod +x "$MM" || { echo "  错误: 无法给 micromamba 加执行权限"; exit 1; }
fi
"$MM" --version | head -1 | sed 's/^/  /'

echo
echo "=== 2. 清华镜像配置 ==="
mkdir -p "$H"
cat > "$H/.condarc" <<'EOF'
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
EOF

echo
echo "=== 3. 创建环境 (python 3.12 + ambertools) ==="
# gmx_MMPBSA 要求 python>=3.11,<3.13, 且需要 cpptraj/tleap/sander(parmed)
# 判断环境是否可用: 以"关键可执行文件都在"为准。
# 注意: 不能只看 `micromamba env list` 的名字 —— 它的输出行首有两个空格,
#       `grep -qE "^mmpbsa"` 匹配不到, 会误判成"环境不存在"从而触发重建。
#       而 micromamba create 会先清空目标环境再下载, 镜像一旦 404 就会留下空壳环境。
ENV_DIR="$H/envs/mmpbsa"
env_ok() {
  [ -x "$ENV_DIR/bin/python" ] && [ -x "$ENV_DIR/bin/gmx_MMPBSA" ] && [ -x "$ENV_DIR/bin/cpptraj" ]
}
if env_ok; then
  echo "  环境 mmpbsa 已存在且完整, 跳过创建"
else
  if [ -d "$ENV_DIR" ]; then
    echo "  环境 mmpbsa 存在但不完整(可能上次创建失败), 先移除"
    "$MM" env remove -n mmpbsa -y >/dev/null 2>&1 || rm -rf "$ENV_DIR"
  fi
  create_env() {
    timeout 2400 "$MM" create -y -n mmpbsa -c conda-forge \
        python=3.12 ambertools=26.0 parmed mpi4py numpy=1.26 pandas matplotlib seaborn scipy tqdm \
        > "$H/create.log" 2>&1
  }
  echo "  创建环境 (python 3.12 + ambertools 26.0) ..."
  if create_env; then
    echo "  create OK"
  else
    echo "  用当前镜像源创建失败, 清理缓存后改用官方 conda-forge 重试 ..."
    tail -5 "$H/create.log" | sed 's/^/    /'
    "$MM" clean --all -y >/dev/null 2>&1 || true
    sed -i 's#https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge#https://conda.anaconda.org/conda-forge#' "$H/.condarc" 2>/dev/null || true
    if create_env; then
      echo "  create OK (官方源)"
    else
      tail -20 "$H/create.log" | sed 's/^/    /'
      echo "  错误: 环境创建失败。日志 $H/create.log"
      exit 1
    fi
  fi
fi

echo
echo "=== 4. 安装 gmx_MMPBSA ==="
if ! env_ok; then
  echo "  错误: mmpbsa 环境不可用, 跳过 gmx_MMPBSA 安装"
  exit 1
fi
"$MM" run -n mmpbsa pip install --no-cache-dir \
    -i "${PIP_INDEX_URL:-https://pypi.tuna.tsinghua.edu.cn/simple}" gmx_MMPBSA 2>&1 | tail -2

echo
echo "=== 5. 验证 ==="
"$MM" run -n mmpbsa bash -c '
  echo "  python : $(python --version 2>&1)"
  for t in cpptraj tleap sander parmchk2 gmx_MMPBSA; do
    printf "  %-12s " "$t"; command -v $t >/dev/null && echo "OK" || echo "缺失"
  done
  python -c "import GMXMMPBSA; print(\"  GMXMMPBSA\", GMXMMPBSA.__version__)"
'

echo
echo "=== 6. GROMACS 需在 PATH 中(给 gmx_MMPBSA 用) ==="
echo "  当前 gmx: $(command -v gmx || echo 未找到)"
echo "  使用时设置: export GMXBIN=/usr/local/bin"
