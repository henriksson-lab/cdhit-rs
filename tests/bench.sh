#!/usr/bin/env bash
# Timing comparison against the C++ reference. Reports the best of N runs of *user CPU
# time*, because the host is shared and wall time is dominated by other tenants.
#
#   BENCH_DNA=<dna.fna> BENCH_AA=<prot.faa> ./tests/bench.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
CPP="${CPP_DIR:-$HERE/../cdhit}"
RS="${RS_DIR:-$HERE/../target/release}"
DNA="${BENCH_DNA:-$HERE/data/dna_mixed.fna}"
AA="${BENCH_AA:-$HERE/data/prot_mixed.faa}"
N="${N:-2}"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

bu() {
  local b=99999 t
  for _ in $(seq "$N"); do
    t=$( { /usr/bin/time -f "%U" "$@" >/dev/null; } 2>&1 | tail -1 )
    b=$(python3 -c "print(min($b,$t))")
  done
  echo "$b"
}

row() {                       # row <label> <prog> <args...>
  local label="$1" prog="$2"; shift 2
  local a b
  a=$(bu "$CPP/$prog" -o "$WORK/c" "$@")
  b=$(bu "$RS/$prog"  -o "$WORK/r" "$@")
  printf '%-26s C++ %7ss   Rust %7ss   %5sx   %s\n' "$label" "$a" "$b" \
    "$(python3 -c "print(f'{$a/$b:.2f}')")" \
    "$(diff -q "$WORK/c.clstr" "$WORK/r.clstr" >/dev/null 2>&1 && echo same || echo DIFF)"
}

echo "user CPU, best of $N; ratio > 1.00 means the Rust port is faster"
row "cd-hit-est -n 5"  cd-hit-est -i "$DNA" -c 0.95 -n 5  -T 1
row "cd-hit-est -n 8"  cd-hit-est -i "$DNA" -c 0.95 -n 8  -T 1
row "cd-hit-est -n 10" cd-hit-est -i "$DNA" -c 0.95 -n 10 -T 1
row "cd-hit-est -n 11" cd-hit-est -i "$DNA" -c 0.95 -n 11 -T 1
row "cd-hit-est c=0.90" cd-hit-est -i "$DNA" -c 0.90 -T 1
row "cd-hit -n 5"      cd-hit     -i "$AA"  -c 0.90 -n 5 -T 1
row "cd-hit -n 4"      cd-hit     -i "$AA"  -c 0.90 -n 4 -T 1
row "cd-hit c=0.70"    cd-hit     -i "$AA"  -c 0.70 -T 1
row "cd-hit-454"       cd-hit-454 -i "$DNA" -c 0.98 -T 1

# 2D takes a second database, so it does not fit `row`.
row2d() {                     # row2d <label> <prog> <db1> <db2> <args...>
  local label="$1" prog="$2" db1="$3" db2="$4"; shift 4
  local a b
  a=$(bu "$CPP/$prog" -i "$db1" -i2 "$db2" -o "$WORK/c" "$@")
  b=$(bu "$RS/$prog"  -i "$db1" -i2 "$db2" -o "$WORK/r" "$@")
  printf '%-26s C++ %7ss   Rust %7ss   %5sx   %s\n' "$label" "$a" "$b" \
    "$(python3 -c "print(f'{$a/$b:.2f}')")" \
    "$(diff -q "$WORK/c.clstr" "$WORK/r.clstr" >/dev/null 2>&1 && echo same || echo DIFF)"
}
row2d "cd-hit-est-2d"  cd-hit-est-2d "$DNA" "$DNA" -c 0.95 -T 1
row2d "cd-hit-2d"      cd-hit-2d     "$AA"  "$AA"  -c 0.90 -T 1
