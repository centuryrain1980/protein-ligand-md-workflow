#!/bin/bash
# skill 回归测试入口
# 用法: bash scripts/tests/run_tests.sh [--fast] [-k 关键字]
set -u
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$DIR/test_all.py" "$@"