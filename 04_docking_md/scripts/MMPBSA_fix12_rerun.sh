#!/bin/bash
# Re-run the three chunks with the patched scripts (no internal PBC clustering).
cd "$(dirname "$0")" || exit 1
for d in fixA fixB fixC; do
  ( cd "$d" && rm -f *_MMPBSA.dat *.pqr *.out && bash gmx_mmpbsa.bash > run.log 2>&1; echo $? > rc.txt ) &
done
wait
echo "ALL DONE"
for d in fixA fixB fixC; do
  echo "$d rc=$(cat $d/rc.txt 2>/dev/null) frames=$(grep -c 'MM-PBSA(with DH)' $d/run.log 2>/dev/null)"
done
