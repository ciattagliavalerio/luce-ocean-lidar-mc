#!/bin/bash
# Usage, from the repository root:
#   runs/2026-11-16_glint_sweep/run_sweep.sh <series_dir> <nphotons> <P>
# Creates one subdirectory per (theta_v, u10) pair, with its input.nml
# derived from input/reference.nml, and runs the simulations P at a time.
# Each subdirectory also gets duration.txt (wall-clock seconds).
set -e
export LC_NUMERIC=C
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"   # Repository root.
OUT="$1"; NPH="$2"; P="$3"
mkdir -p "$OUT"
for a in 0 6 12 20 30; do
  for u in 0 2 5 7.5 10 15; do
    d=$(printf "a%02d_U%04.1f" "$a" "$u")
    mkdir -p "$OUT/$d"
    sed -E \
      -e "s/^(  theta_v +=) 12\.0 /\1 $(printf '%-4s' "$a.0") /" \
      -e "s/^(  u10 +=) 7\.5 /\1 $(printf '%-3s' "$u") /" \
      -e "s/^(  nphotons =) 10000000 /\1 $(printf '%-8s' "$NPH") /" \
      "$ROOT/input/reference.nml" > "$OUT/$d/input.nml"
  done
done
cd "$OUT"
ls -d a*_U*/ | xargs -P "$P" -I{} bash -c \
  'cd {} && s=$(date +%s.%N) && "'"$ROOT"'/mc_water" input.nml > run.log 2>&1; e=$(date +%s.%N); echo "$e - $s" | bc > duration.txt'
