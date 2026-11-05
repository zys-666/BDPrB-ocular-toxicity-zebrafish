#!/bin/bash
cd "$(dirname "$0")" || exit 1
tail -4 driver.log
echo "--- frames done ---"
for d in fixA fixB fixC; do
  n=$(grep -c 'MM-PBSA(with DH)' "$d/run.log" 2>/dev/null)
  echo "$d: ${n:-0}"
done
