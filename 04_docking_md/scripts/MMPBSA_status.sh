#!/bin/bash
cd /mnt/d/Gromacs25_1_win/mmpbsa/gmxmmpbsa
echo "=== 当前结果文件 ==="
for f in _pid~MMPBSA.dat _pid~res.dat _pid~res_MMPBSA.dat _pid~resMM.dat _pid~resMM_COU.dat _pid~resMM_VDW.dat _pid~resPBSA.dat _pid~resPBSA_PB.dat _pid~resPBSA_SA.dat; do
  if [ -f "$f" ]; then
    printf '  %-24s %8d B\n' "$f" "$(stat -c%s "$f")"
  else
    printf '  %-24s <缺失>\n' "$f"
  fi
done
echo
echo "=== _pid~MMPBSA.dat 内容判定 ==="
echo "  样本: $(head -2 _pid~MMPBSA.dat | tail -1 | cut -c1-80)"
echo "  帧数: $(grep -cE '^_pid~[0-9.]+ns' _pid~MMPBSA.dat)"
echo "  首帧时间: $(awk '/^_pid~/{print $1; exit}' _pid~MMPBSA.dat)"
echo "  末帧时间: $(grep -oE '_pid~[0-9.]+ns' _pid~MMPBSA.dat | tail -1)"
echo
echo "=== 其他可用备份 ==="
ls -la old_run_backup/ 2>/dev/null | tail -n +4
