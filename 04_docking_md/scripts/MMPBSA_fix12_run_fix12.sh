#!/bin/bash
# Re-run MM-PBSA for the 12 frames that failed PBC handling in the segmented run.
set -e
cd "$(dirname "$0")"

GMX=/mnt/d/Gromacs25_1_win/gromacs/bin/gmx.exe
SRC=../../MD/03_production/md_whole.xtc
TPR=../md.tpr
NDX=../index.ndx
BASE=../gmx_mmpbsa.bash

mkdir -p frames
TIMES="58400 58700 59300 60500 60800 87000 87600 87900 88200 89700 95300 96500"
rm -f frames/f_*.xtc

for t in $TIMES; do
  printf 'Protein_LIGAND\n' | "$GMX" trjconv -f "$SRC" -s "$TPR" -n "$NDX" \
      -dump "$t" -o "frames/f_${t}.xtc" >/dev/null 2>&1
done
ls frames/f_*.xtc | sort -V > frame_list.txt
wc -l < frame_list.txt

"$GMX" trjcat -f $(cat frame_list.txt | tr '\n' ' ') -o fix12.xtc -cat >/dev/null 2>&1
printf 'Protein_LIGAND\n' | "$GMX" check -f fix12.xtc 2>&1 | grep -E 'Step|Last frame' | head -3

# three chunks (by original time)
printf 'Protein_LIGAND\n' | "$GMX" trjconv -f fix12.xtc -s "$TPR" -n "$NDX" -b 58000 -e 61000 -o chunkA.xtc >/dev/null 2>&1
printf 'Protein_LIGAND\n' | "$GMX" trjconv -f fix12.xtc -s "$TPR" -n "$NDX" -b 86000 -e 89000 -o chunkB.xtc >/dev/null 2>&1
printf 'Protein_LIGAND\n' | "$GMX" trjconv -f fix12.xtc -s "$TPR" -n "$NDX" -b 89000 -e 97000 -o chunkC.xtc >/dev/null 2>&1

for tag in A B C; do
  d="fix$tag"
  rm -rf "$d"; mkdir -p "$d"
  cp "$BASE" "$d/gmx_mmpbsa.bash"
  cp "$TPR" "$d/md.tpr"
  cp "$NDX" "$d/index.ndx"
  cp "chunk$tag.xtc" "$d/chunk$tag.xtc"
  sed -i "s|^trj=.*|trj=chunk$tag.xtc|" "$d/gmx_mmpbsa.bash"
  sed -i "s|^pid=.*|pid=fix$tag|" "$d/gmx_mmpbsa.bash"
  # remove the hard-coded -dt 300 so irregular frames are kept
  sed -i 's|^trjconv=.*|trjconv="$gmx trjconv"|' "$d/gmx_mmpbsa.bash"
  ( cd "$d" && bash gmx_mmpbsa.bash > run.log 2>&1; echo $? > rc.txt ) &
done
wait
echo "ALL DONE"
for tag in A B C; do
  echo "fix$tag rc=$(cat fix$tag/rc.txt 2>/dev/null) frames=$(grep -c 'MM-PBSA(with DH)' fix$tag/run.log 2>/dev/null)"
done
