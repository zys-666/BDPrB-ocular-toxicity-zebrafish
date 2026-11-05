#!/bin/bash
################################################################################
#  MMPBSA 并行驱动程序  (帧数无关 / frame-count agnostic)
#
#  原理:
#    gmx_mmpbsa 每帧独立计算, 帧间无依赖 -> 可并行
#    本脚本用 gmx trjconv 的 -b/-e 把轨迹按时间切成 N 段,
#    每段起一个独立的 gmx_mmpbsa 进程(独立的 pid 前缀, 输出不冲突),
#    全部跑完后自动合并成一份完整的 _pid~*.dat
#
#  关键设计 (为了"帧数无关"):
#    1. 不按"帧数"切, 按 trjconv 的 -b/-e 时间边界切 -> 帧数任意都行
#    2. 启动前先探测轨迹的真实时间范围与步长, 自动均分
#    3. 合并时按时间戳排序, 不依赖帧数
#
#  用法:
#    bash run_mmpbsa_parallel.sh                       # 全自动: 探测+均分+跑+合并
#    bash run_mmpbsa_parallel.sh -j 8                  # 指定 8 路并行
#    bash run_mmpbsa_parallel.sh -f traj.xtc -s a.tpr -n index.ndx -pro Protein -lig MOL
#    bash run_mmpbsa_parallel.sh -t "50000 100000"     # 指定时间窗口 (ps)
#    bash run_mmpbsa_parallel.sh -dry                  # 只显示切分方案, 不跑
#
#  依赖: 已修好的 gmx_mmpbsa.bash (WSL 64位 APBS 版)
################################################################################
set -u

############ 可调默认值 ############
WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MMPBSA_SCRIPT="$WORKDIR/gmx_mmpbsa.bash"     # 主脚本
# 兼容被重命名的版本
[ -f "$MMPBSA_SCRIPT" ] || MMPBSA_SCRIPT="$WORKDIR/gmx_mmpbsa_WSL.bash"
JOBS=auto                                     # 并行路数 (auto = 按内存/核数自动)
TRJ="trj_70-100ns.xtc"
TPR="md_0_1.tpr"
NDX="index.ndx"
PRO="Protein"
LIG="	MOL"
DT=300                                        # 采样步长 ps (0=用轨迹原始步长)
TIME_WINDOW=""                                # "begin end" ps, 空=用整条轨迹
DRY=0
MEM_PER_JOB_GB=3                              # 单个 APBS 进程预估峰值内存
KEEP_TMP=0                                    # 1=保留分段中间文件

############ 参数解析 ############
while [ $# -gt 0 ]; do
  case "$1" in
    -j)   JOBS="$2"; shift 2;;
    -f)   TRJ="$2"; shift 2;;
    -s)   TPR="$2"; shift 2;;
    -n)   NDX="$2"; shift 2;;
    -pro) PRO="$2"; shift 2;;
    -lig) LIG="$2"; shift 2;;
    -dt)  DT="$2"; shift 2;;
    -t)   TIME_WINDOW="$2"; shift 2;;
    -mem) MEM_PER_JOB_GB="$2"; shift 2;;
    -dry) DRY=1; shift;;
    -keep) KEEP_TMP=1; shift;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) echo "未知参数: $1"; exit 1;;
  esac
done

cd "$WORKDIR" || { echo "!! 无法进入 $WORKDIR"; exit 1; }

############ 环境 ############
source /root/miniforge3/etc/profile.d/conda.sh
conda activate apbs64
GMX=/mnt/d/Gromacs25_1_win/gromacs/bin/gmx.exe
export MCSH_HOME="$WORKDIR"

echo "=============================================================="
echo " MMPBSA 并行驱动"
echo "=============================================================="
echo " 工作目录 : $WORKDIR"
echo " 主脚本   : $MMPBSA_SCRIPT"
echo " 轨迹     : $TRJ"
echo " 拓扑     : $TPR"
echo " 索引     : $NDX"
echo " 进程组   : $PRO / $LIG"

# 检查文件
for f in "$MMPBSA_SCRIPT" "$TRJ" "$TPR" "$NDX"; do
  [ -f "$f" ] || { echo "!! 缺文件: $f"; exit 1; }
done

############ 1. 探测轨迹真实时间范围 ############
echo
echo "[1] 探测轨迹时间范围..."
# gmx check 的信息走 stderr, 需 2>&1
# 组合组名 = Protein_MOL (脚本内部会自己建, 这里只用默认组探测)
PROBE=$(printf '0\n' | "$GMX" check -f "$TRJ" 2>&1)

# 帧数与步长: 找 "Step   <N>   <dt>" 行
NFRM=$(echo "$PROBE"  | awk '/^Step/{print $2; exit}')
TSTEP=$(echo "$PROBE" | awk '/^Step/{print $3; exit}')
# 首末帧时间: 找 "Reading frame   N time  T"  (注意 CRLF 和多余空格)
FIRST_T=$(echo "$PROBE" | tr -d '\r' | grep 'Reading frame' | head -1 | awk '{print $NF}')
LAST_T=$(echo  "$PROBE" | tr -d '\r' | grep 'Reading frame' | tail -1 | awk '{print $NF}')

# 若 Reading frame 被省略(大轨迹, gmx 只每10帧打印), 用帧数*步长估算
if [ -z "${FIRST_T:-}" ] || ! echo "$FIRST_T" | grep -qE '^[0-9.]+$'; then
  FIRST_T=$(awk -v d="$TSTEP" 'BEGIN{print d}')
fi
if [ -z "${LAST_T:-}" ] || ! echo "$LAST_T" | grep -qE '^[0-9.]+$'; then
  LAST_T=$(awk -v n="$NFRM" -v d="$TSTEP" 'BEGIN{print n*d}')
fi

if [ -z "${NFRM:-}" ] || [ -z "${TSTEP:-}" ]; then
  echo "!! 无法探测轨迹信息, 请手动指定 -t \"begin end\""
  exit 1
fi

echo "   帧数     : $NFRM"
echo "   时间步长 : $TSTEP ps"
echo "   首帧时间 : ${FIRST_T:-?} ps"
echo "   末帧时间 : ${LAST_T:-?} ps"

# 确定时间窗口
if [ -n "$TIME_WINDOW" ]; then
  TB=$(echo "$TIME_WINDOW" | awk '{print $1}')
  TE=$(echo "$TIME_WINDOW" | awk '{print $2}')
else
  TB=${FIRST_T:-0}
  TE=${LAST_T:-$(awk -v n="$NFRM" -v d="$TSTEP" 'BEGIN{print (n-1)*d}')}
fi
# 去掉小数
TB=$(printf '%.0f' "$TB"); TE=$(printf '%.0f' "$TE")
SPAN=$((TE - TB))
echo "   计算窗口 : $TB -> $TE ps  (跨度 ${SPAN} ps)"

# 采样步长
[ "$DT" -le 0 ] && DT="$TSTEP"
N_EFF=$(( SPAN / DT + 1 ))
echo "   采样步长 : $DT ps  ->  约 $N_EFF 帧参与计算"

############ 2. 决定并行路数 ############
if [ "$JOBS" = "auto" ]; then
  NCPU=$(nproc)
  MEM_AVAIL_GB=$(free -g | awk 'NR==2{print $7}')
  BY_MEM=$(awk -v m="$MEM_AVAIL_GB" -v p="$MEM_PER_JOB_GB" 'BEGIN{j=int((m-4)/p); print (j<1?1:j)}')  # 留 4G 给系统
  JOBS=$BY_MEM
  [ "$JOBS" -gt "$NCPU" ] && JOBS=$NCPU
  [ "$JOBS" -gt "$N_EFF" ] && JOBS=$N_EFF
  echo
  echo "[2] 自动决定并行路数"
  echo "    CPU 核数      : $NCPU"
  echo "    可用内存      : ${MEM_AVAIL_GB} GB  (每进程按 ${MEM_PER_JOB_GB} GB 估)"
  echo "    内存允许      : $BY_MEM 路"
  echo "    -> 采用       : $JOBS 路"
else
  echo
  echo "[2] 手动指定并行路数: $JOBS"
  [ "$JOBS" -gt "$N_EFF" ] && { echo "    (帧数只有 $N_EFF, 自动降到 $N_EFF)"; JOBS=$N_EFF; }
fi
[ "$JOBS" -lt 1 ] && JOBS=1

############ 3. 计算时间切分点 ############
# 把 [TB, TE] 均分成 JOBS 段, 每段用 -b/-e 指定
echo
echo "[3] 时间切分方案"
SEG_B=(); SEG_E=()
STEP=$(( SPAN / JOBS ))
for ((i=0; i<JOBS; i++)); do
  b=$(( TB + i*STEP ))
  if [ $i -eq $((JOBS-1)) ]; then e=$TE; else e=$(( TB + (i+1)*STEP )); fi
  SEG_B+=("$b"); SEG_E+=("$e")
  n=$(( (e-b)/DT + 1 ))
  printf "    段%-2d : %8d -> %8d ps   约 %4d 帧\n" "$((i+1))" "$b" "$e" "$n"
done

if [ "$DRY" -eq 1 ]; then
  echo
  echo "[dry-run] 不执行, 仅显示方案"
  exit 0
fi

############ 4. 并行执行 ############
echo
echo "[4] 启动 $JOBS 路并行..."
echo "=============================================================="

RUN_NAME="par$(date +%m%d_%H%M%S)"
mkdir -p "$RUN_NAME"
LOGDIR="$WORKDIR/$RUN_NAME"
PIDS=()

for ((i=0; i<JOBS; i++)); do
  seg=$((i+1))
  b=${SEG_B[$i]}; e=${SEG_E[$i]}

  # 每段独立的 pid 前缀 -> 输出文件天然隔离
  # 用 trjconv 切出该时间段的子轨迹 (组名用 ${PRO}_${LIG}, 与主脚本一致)
  SEGGRP="${PRO}_${LIG}"
  if ! grep -q "\[ *${SEGGRP} *\]" "$NDX" 2>/dev/null; then
    # 索引里没有组合组, 退回第一个可用组
    SEGGRP=$(grep -m1 -oP '(?<=\[ ).*(?= \])' "$NDX" | head -1)
    [ -z "$SEGGRP" ] && SEGGRP="System"
  fi
  seg_trj="seg${seg}_${TRJ}"
  printf '%s\n' "$SEGGRP" | "$GMX" trjconv -s "$TPR" -n "$NDX" \
      -f "$TRJ" -o "$seg_trj" -b "$b" -e "$e" -dt "$DT" >/dev/null 2>&1
  if [ ! -f "$seg_trj" ]; then
    echo "    !! 段$seg 切分失败 ($SEGGRP, $b-$e ps)"
    continue
  fi

  sub="$LOGDIR/seg${seg}"
  mkdir -p "$sub"
  # 注意: Windows 的 gmx.exe 无法跟随 WSL 符号链接(error 22),
  #       所以这里必须【复制】而不是 ln -s
  cp "$WORKDIR/$seg_trj" "$sub/$seg_trj"
  cp "$WORKDIR/$TPR"     "$sub/$TPR"
  cp "$WORKDIR/$NDX"     "$sub/$NDX"
  cp "$MMPBSA_SCRIPT"    "$sub/gmx_mmpbsa.bash"

  (
    cd "$sub"
    # 把脚本里的 trj 改成这一段
    sed -i "s|^trj=.*|trj=$seg_trj|" gmx_mmpbsa.bash
    sed -i "s|^pid=.*|pid=seg${seg}|" gmx_mmpbsa.bash
    bash gmx_mmpbsa.bash > "seg${seg}.log" 2>&1
    echo "$?" > "seg${seg}.rc"
  ) &
  PIDS+=($!)
  printf "    段%-2d (%d-%d ps) 已启动, PID=%d\n" "$seg" "$b" "$e" "$!"
done

echo
echo "    等待全部完成 (Ctrl-C 可中断, 已完成的段会保留)..."
FAIL=0
for i in "${!PIDS[@]}"; do
  seg=$((i+1))
  if wait "${PIDS[$i]}"; then
    rc=$(cat "$LOGDIR/seg${seg}/seg${seg}.rc" 2>/dev/null || echo 0)
    if [ "$rc" = "0" ]; then
      echo "    ✅ 段$seg 完成"
    else
      echo "    ⚠️  段$seg 退出码=$rc (检查 $LOGDIR/seg${seg}/seg${seg}.log)"
      FAIL=$((FAIL+1))
    fi
  else
    echo "    ❌ 段$seg 异常"
    FAIL=$((FAIL+1))
  fi
done

echo "=============================================================="
echo "[5] 合并结果..."
bash "$WORKDIR/merge_parallel.sh" "$LOGDIR" "$WORKDIR"
MRC=$?

# 若因已存在结果文件而中止, 自动改用带时间戳的安全目录, 不覆盖
if [ "$MRC" -ne 0 ]; then
  echo
  echo "  (合并未完成; 若因目标已存在被拒, 下面自动改用新目录重试一次)"
  ALT="$LOGDIR/merged"
  mkdir -p "$ALT"
  if bash "$WORKDIR/merge_parallel.sh" "$LOGDIR" "$ALT"; then
    echo "  ✅ 已合并到: $ALT"
    echo "     (原目录结果未被改动)"
    MRC=0
  fi
fi

echo
echo "=============================================================="
if [ "$FAIL" -eq 0 ] && [ "$MRC" -eq 0 ]; then
  echo " ✅ 全部完成"
else
  echo " ⚠️  有 $FAIL 段失败, 合并返回 $MRC"
fi
echo " 分段目录: $LOGDIR"
echo " 合并结果: $WORKDIR/_pid~*.dat"
echo "=============================================================="
