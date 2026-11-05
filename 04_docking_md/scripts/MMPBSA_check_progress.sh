#!/bin/bash
cd "$(dirname "$0")" || exit 1
d=$(ls -d par* | tail -1)
echo "DIR=$d"
tot=0
for f in "$d"/seg*/seg*.log; do
  n=$(grep -c 'MM-PBSA(with DH)' "$f" 2>/dev/null)
  s=$(echo "$f" | sed 's|.*/\(seg[0-9]*\)/.*|\1|')
  echo "  $s : $n/28 frames"
  tot=$((tot+n))
done
echo "TOTAL: $tot / 167"
