#!/bin/bash
# 安装 gmx_MMPBSA 运行环境(Python 3.12 + AmberTools + gmx_MMPBSA)
# 用法: bash setup_mmpbsa.sh
set -u
H="${MMPBSA_HOME:-$HOME/micromamba}"
export MAMBA_ROOT_PREFIX="$H"
MM="$H/bin/micromamba"

echo "=== 1. micromamba ==="
if [ ! -x "$MM" ]; then
  mkdir -p "$H"
  echo "  请把 micromamba.tar.bz2 放到工作区, 或从 conda-forge 下载:"
  echo "    https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/micromamba-2.8.1-1.tar.bz2"
  TAR=$(ls -1 "$HOME"/../*/work/ref/micromamba.tar.bz2 2>/dev/null | head -1)
  [ -z "$TAR" ] && TAR=$(find /mnt -maxdepth 6 -name "micromamba.tar.bz2" 2>/dev/null | head -1)
  [ -z "$TAR" ] && { echo "  错误: 找不到 micromamba.tar.bz2"; exit 1; }
  echo "  从 $TAR 解压"
  tar xjf "$TAR" -C "$H" bin/micromamba
  chmod +x "$MM"
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
if "$MM" env list 2>/dev/null | grep -qE "^mmpbsa[[:space:]]"; then
  echo "  环境 mmpbsa 已存在, 跳过创建"
else
  timeout 2400 "$MM" create -y -n mmpbsa -c conda-forge \
      python=3.12 ambertools=26.0 parmed mpi4py numpy=1.26 pandas matplotlib seaborn scipy tqdm \
      > "$H/create.log" 2>&1
  echo "  create exit=$?  (日志 $H/create.log)"
fi

echo
echo "=== 4. 安装 gmx_MMPBSA ==="
"$MM" run -n mmpbsa pip install --no-cache-dir \
    -i https://pypi.tuna.tsinghua.edu.cn/simple gmx_MMPBSA 2>&1 | tail -2

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
