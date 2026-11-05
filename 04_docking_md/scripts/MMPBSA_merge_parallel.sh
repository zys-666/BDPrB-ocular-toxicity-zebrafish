#!/bin/bash
################################################################################
#  合并 N 路并行 MMPBSA 的结果
#
#  用法: bash merge_parallel.sh <分段目录> <输出目录>
#
#  原理:
#    - 每段输出 seg<N>~*.dat, 内含该段的逐帧/逐残基数据
#    - 逐帧文件: 按帧时间戳排序后串接
#    - 逐残基文件: 数值需按帧数加权平均 (每段帧数不同!)
#    - 尾部摘要: 用合并后的全量数据重新计算
#
#  帧数无关: 全部靠解析时间戳, 不依赖任何硬编码帧数
################################################################################
set -u

LOGDIR="${1:?用法: merge_parallel.sh <分段目录> <输出目录>}"
OUTDIR="${2:?缺少输出目录}"

# 先用绝对路径固化两个目录 (因为下面要 cd, 相对路径会失效)
SCRIPT_CWD=$(pwd)
LOGDIR=$(cd "$LOGDIR" 2>/dev/null && pwd) || { echo "!! 分段目录不存在: $1"; exit 1; }
mkdir -p "$OUTDIR" || { echo "!! 无法创建输出目录: $OUTDIR"; exit 1; }
OUTDIR=$(cd "$OUTDIR" && pwd)

cd "$OUTDIR" || exit 1

echo "  分段目录: $LOGDIR"
echo "  输出目录: $OUTDIR"

# ===== 安全保护: 拒绝静默覆盖已有结果 =====
if [ "${FORCE:-0}" != "1" ]; then
  CLASH=""
  for f in _pid~MMPBSA.dat _pid~res.dat; do
    [ -f "$f" ] && CLASH="$CLASH $f"
  done
  if [ -n "$CLASH" ]; then
    echo
    echo "  !! 输出目录已存在结果文件:$CLASH"
    echo "     为避免覆盖, 已中止。"
    echo "     如确认要覆盖:  FORCE=1 bash merge_parallel.sh $LOGDIR $OUTDIR"
    echo "     或改用其他输出目录: bash merge_parallel.sh $LOGDIR /path/to/newdir"
    exit 1
  fi
fi
# ===== 保护结束 =====

# 找出所有段 (LOGDIR 已是绝对路径)
SEGS=$(ls -d "$LOGDIR"/seg* 2>/dev/null | sort -V)
NSEG=$(echo "$SEGS" | grep -c . )
echo "  发现 $NSEG 个段"

if [ "$NSEG" -lt 1 ]; then echo "  !! 没有找到任何分段结果"; exit 1; fi

############ 工具: 从逐帧 .dat 提取 (时间戳, 整行) ############
# 逐帧行形如:  seg1~80.0ns   -28.249(  -26.044) | ...
# 要提取出 "80.0" 作为排序键
extract_keyed() {
  local f="$1"
  # 帧行前缀可能是 seg1~ / _seg1~ / _pid~ / par1~ / pid~
  # 统一策略: 抓时间键, 并把首字段【规范化】为 _pid~<t>ns, 使合并结果与单进程版完全一致
  # 同时把表头行(#Frame ...)单独输出到 header 通道, 供合并时使用
  awk '
    # 表头: 含 "Binding" 或首个非空字段以 # 开头 -> 写到 header 文件
    /Binding|^ *#Frame/ && !hdr_done {
      print > HEADER_FILE
      hdr_done = 1
      next
    }
    {
      t = $1
      sub(/.*~/, "", t)      # 去前缀 -> 80.0ns
      sub(/ns$/, "", t)      # -> 80.0
      if (t ~ /^[0-9]+(\.[0-9]+)?$/) {
        $1 = "_pid~" t "ns"   # 规范化前缀
        printf "%012.4f\t%s\n", t+0, $0
      }
    }
  ' HEADER_FILE="$2" "$f"
}

############ 1. 合并逐帧总表 _pid~MMPBSA.dat ############
echo
echo "  [1] 合并逐帧总表..."

TMP=$(mktemp); HDRF=$(mktemp)
# 表头由 extract_keyed 从分段文件里原样抓取, 写入 HDRF
for s in $SEGS; do
  for f in "$s"/*~MMPBSA.dat; do
    [ -f "$f" ] || continue
    extract_keyed "$f" "$HDRF" >> "$TMP"
  done
done
# 按时间排序, 去重 (不同段边界可能重叠), 还原成原格式
{
  # 表头 (原样保留, 含前导空格)
  [ -s "$HDRF" ] && head -1 "$HDRF"
  sort -t$'\t' -k1,1n "$TMP" | awk -F'\t' '!seen[$1]++ {print $2}'
} > "$OUTDIR/_pid~MMPBSA.dat"
rm -f "$HDRF"
# 统计合并后的帧数 (任意前缀)
NF=$(awk '{t=$1; sub(/.*~/,"",t); sub(/ns$/,"",t); if(t ~ /^[0-9]+(\.[0-9]+)?$/) n++} END{print n+0}' "$OUTDIR/_pid~MMPBSA.dat")
echo "      逐帧总表: $NF 帧"

############ 2. 合并逐帧逐残基表 ############
echo "  [2] 合并逐帧逐残基表..."
for suffix in res_MMPBSA resMM resMM_VDW resMM_COU resPBSA resPBSA_PB resPBSA_SA; do
  TMP=$(mktemp); HDR=""; HDRF=$(mktemp)
  found=0
  for s in $SEGS; do
    for f in "$s"/*~${suffix}.dat; do
      [ -f "$f" ] || continue
      found=1
      extract_keyed "$f" "$HDRF" >> "$TMP"
    done
  done
  # 取抓到的第一个表头 (原样, 含前导空格)
  [ -s "$HDRF" ] && HDR=$(head -1 "$HDRF")
  rm -f "$HDRF"
  if [ "$found" -eq 1 ] && [ -s "$TMP" ]; then
    {
      [ -n "$HDR" ] && echo "$HDR"
      sort -t$'\t' -k1,1n "$TMP" | awk -F'\t' '!seen[$1]++ {print $2}'
    } > "$OUTDIR/_pid~${suffix}.dat"
    n=$(awk '{t=$1; sub(/.*~/,"",t); sub(/ns$/,"",t); if(t ~ /^[0-9]+(\.[0-9]+)?$/) c++} END{print c+0}' "$OUTDIR/_pid~${suffix}.dat")
    echo "      _pid~${suffix}.dat : $n 帧"
  fi
  rm -f "$TMP"
done

############ 3. 重新计算逐残基平均表 _pid~res.dat ############
# 关键: 不能简单拼各段的平均, 必须按各段帧数加权
echo "  [3] 重算逐残基平均表 (按帧数加权)..."

python3 - "$LOGDIR" "$OUTDIR" <<'PYEOF'
import sys, os, re, glob, math
logdir, outdir = sys.argv[1], sys.argv[2]

segs = sorted(glob.glob(os.path.join(logdir, 'seg*')), key=lambda p: int(re.search(r'seg(\d+)$', p).group(1)))
if not segs:
    print("      !! 无分段"); sys.exit(1)

# 每段的帧数 = 该段 *_MMPBSA.dat 里符合 "<前缀>~<数字>ns" 的数据行数
FRAME_RE = re.compile(r'^\S*~\d+(\.\d+)?ns\b')
seg_nframes = {}
for s in segs:
    n = 0
    for f in glob.glob(os.path.join(s, '*~MMPBSA.dat')):
        for l in open(f, errors='replace'):
            if FRAME_RE.match(l.strip()):
                n += 1
    seg_nframes[s] = n
total = sum(seg_nframes.values())
print("      各段帧数: %s   合计 %d" % (
    {os.path.basename(k): v for k, v in seg_nframes.items()}, total))
if total == 0:
    print("      !! 总帧数为 0"); sys.exit(1)

# 累加各段的 res.dat (格式: Lig/Pro <res> + 10 个数值)
acc = {}          # tag -> [sum_i]
order = []        # 保持残基顺序
hdr = None

for s in segs:
    for f in glob.glob(os.path.join(s, '*~res.dat')):
        lines = [l.rstrip('\n') for l in open(f, errors='replace')]
        if hdr is None and lines and lines[0].lstrip().startswith('#'):
            hdr = lines[0]
        for l in lines:
            m = re.match(r'^(Lig|Pro)\s+(\S+)\s+(.*)$', l.strip())
            if not m: continue
            tag = "%s %s" % (m.group(1), m.group(2))
            vals = [float(x) for x in re.findall(r'-?\d+\.\d+', m.group(3))]
            if tag not in acc:
                acc[tag] = [0.0]*len(vals)
                order.append(tag)
            for i, v in enumerate(vals):
                acc[tag][i] += v

# 加权平均 (每段的 res.dat 是该段帧数的平均, 权重 = 该段帧数)
# 但各段 res.dat 已经各自平均, 且各段帧数已知 -> 加权
wmap = {}
for s in segs:
    for f in glob.glob(os.path.join(s, '*~res.dat')):
        # 该段的帧数
        wmap[os.path.dirname(f)] = seg_nframes.get(os.path.dirname(f), 0)

# 重新按 (段帧数, 段值) 加权
weighted = {}
for tag in acc:
    weighted[tag] = None
for s in segs:
    w = seg_nframes.get(s, 0)
    if w == 0: continue
    for f in glob.glob(os.path.join(s, '*~res.dat')):
        for l in open(f, errors='replace'):
            m = re.match(r'^(Lig|Pro)\s+(\S+)\s+(.*)$', l.strip())
            if not m: continue
            tag = "%s %s" % (m.group(1), m.group(2))
            vals = [float(x) for x in re.findall(r'-?\d+\.\d+', m.group(3))]
            if weighted[tag] is None:
                weighted[tag] = [v*w for v in vals]
            else:
                for i, v in enumerate(vals):
                    weighted[tag][i] += v*w

# 输出
def sub_of(tag):
    return '1MOL' if tag.startswith('Lig') else tag.split()[1]

out = os.path.join(outdir, '_pid~res.dat')
fmt = "%-12s %9.3f(%9.3f) %9.3f(%9.3f) %9.3f %9.3f(%9.3f) %9.3f %9.3f %9.3f\n"
with open(out, 'w') as w:
    if hdr: w.write(hdr + "\n")
    for tag in order:
        v = weighted[tag]
        if v is None: continue
        v = [x/total for x in v]
        label = ("%-4s %s" % (tag.split()[0], tag.split()[1]))
        w.write(fmt % (label, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9]))
print("      写出 _pid~res.dat : %d 行 (%d 残基 + 配体)" % (len(order), len(order)-1))
PYEOF

############ 4. 重算尾部摘要 ############
echo "  [4] 重算汇总摘要 (dH / -TdS / dG / Ki)..."
python3 - "$OUTDIR" <<'PYEOF'
import sys, os, re, math
outdir = sys.argv[1]
f = os.path.join(outdir, '_pid~MMPBSA.dat')
FRAME_RE = re.compile(r'^\S*~\d+(\.\d+)?ns\b')
rows = []
for l in open(f, errors='replace'):
    if FRAME_RE.match(l.strip()):
        segs = [t.strip() for t in l.split('ns',1)[-1].split('|') if t.strip()]
        v = []
        for t in segs:
            v += [float(x) for x in re.findall(r'-?\d+\.\d+', t)]
        if len(v) >= 15:
            rows.append(v)
if not rows:
    print("      !! 无数据"); sys.exit(1)

N = len(rows)
import statistics as st
def col(i): return [r[i] for r in rows]

dH  = sum(col(0))/N
MM  = sum(col(2))/N
PB  = sum(col(4))/N
SA  = sum(col(5))/N
COU = sum(col(6))/N
VDW = sum(col(8))/N

RT2kJ = 8.314462618*298.15/1E3
# 相互作用熵 (IE): 用 MM 列
TdS = -RT2kJ*math.log(sum(math.exp((r[2]-MM)/RT2kJ) for r in rows)/N)
dG  = dH - TdS
Ki  = math.exp(dG/RT2kJ)

with open(f, 'a') as w:
    w.write("----------------------------------|------------------------------------------|--------------------------------|\n")
    w.write("%-12s %9.3f(%9.3f) | %9.3f(%9.3f) %9.3f %9.3f | %9.3f(%9.3f) %9.3f |\n" % (
        "  dH=", dH, dH, MM, MM, PB, SA, COU, COU, VDW))
    w.write("%-12s %9.3f(%9.3f) |\n" % ("-TdS=", -TdS, -TdS))
    w.write("----------------------------------|\n")
    w.write("%-12s %9.3f(%9.3f) kJ/mol = %9.3f(%9.3f) kcal/mol\n" % (
        "  dG=", dG, dG, dG/4.184, dG/4.184))
    w.write("%-12s %9.3E(%9.3E)  uM    = %9.3E(%9.3E)    nM\n" % (
        "  Ki=", Ki*1E6, Ki*1E6, Ki*1E9, Ki*1E9))
    w.write("+==============================================================+\n")

print("      dH = %.3f kJ/mol  (%.3f kcal/mol)" % (dH, dH/4.184))
print("      dG = %.3f kJ/mol  (%.3f kcal/mol)  [含 IE 熵]" % (dG, dG/4.184))
print("      合并帧数 = %d" % N)
PYEOF

echo "  ✅ 合并完成"
