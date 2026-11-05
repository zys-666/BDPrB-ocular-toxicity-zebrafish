#!/bin/bash
cd "$(dirname "$0")" || exit 1
tail -6 run_mmpbsa.log
echo "--- segments ---"
d=$(ls -d par* 2>/dev/null | tail -1)
echo "dir=$d"
for s in "$d"/seg*; do
  n=$(grep -c 'MM-PBSA(with DH)' "$s"/*.log 2>/dev/null)
  echo "$s : ${n:-0} frames"
done
