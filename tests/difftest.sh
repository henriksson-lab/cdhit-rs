#!/usr/bin/env bash
# Differential test: run the C++ cd-hit / cd-hit-est and the Rust port over the same
# inputs and options, and require byte-identical .clstr and FASTA output.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
CPP="${CPP_DIR:-$HERE/../cdhit}"
RS="${RS_DIR:-$HERE/../target/release}"
DATA="$HERE/data"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

pass=0; fail=0; failed=()

run_case() {                      # run_case <prog> <input> <label> <opts...>
  local prog="$1" input="$2" label="$3"; shift 3
  local c="$WORK/c.$label" r="$WORK/r.$label"
  "$CPP/$prog" -i "$input" -o "$c" "$@" > "$WORK/c.log" 2>&1
  local crc=$?
  "$RS/$prog"  -i "$input" -o "$r" "$@" > "$WORK/r.log" 2>&1
  local rrc=$?
  if [ $crc -ne $rrc ]; then
    echo "FAIL $label  (exit $crc vs $rrc)"; sed -n '1,15p' "$WORK/r.log"
    fail=$((fail+1)); failed+=("$label"); return
  fi
  if [ ! -f "$c.clstr" ] && [ ! -f "$r.clstr" ]; then
    # both refused the run; require the same diagnostic
    if diff -q <(grep -a "Fatal Error" -A2 "$WORK/c.log") \
               <(grep -a "Fatal Error" -A2 "$WORK/r.log") >/dev/null 2>&1; then
      pass=$((pass+1))
    else
      echo "FAIL $label  (different refusal)"
      diff <(grep -a "Fatal Error" -A2 "$WORK/c.log") <(grep -a "Fatal Error" -A2 "$WORK/r.log") | head -8
      fail=$((fail+1)); failed+=("$label")
    fi
    return
  fi
  if diff -q "$c.clstr" "$r.clstr" >/dev/null 2>&1 && diff -q "$c" "$r" >/dev/null 2>&1; then
    pass=$((pass+1))
  else
    echo "FAIL $label"
    diff "$c.clstr" "$r.clstr" 2>&1 | head -8
    diff -q "$c" "$r" 2>&1 | head -2
    fail=$((fail+1)); failed+=("$label")
  fi
}

est()  { run_case cd-hit-est "$@"; }
prot() { run_case cd-hit     "$@"; }
p454() { run_case cd-hit-454 "$@"; }

# 2D: db2 is compared against db1, which is never clustered.
run_2d() {                        # run_2d <prog> <db1> <db2> <label> <opts...>
  local prog="$1" db1="$2" db2="$3" label="$4"; shift 4
  local c="$WORK/c.$label" r="$WORK/r.$label"
  "$CPP/$prog" -i "$db1" -i2 "$db2" -o "$c" "$@" > "$WORK/c.log" 2>&1
  "$RS/$prog"  -i "$db1" -i2 "$db2" -o "$r" "$@" > "$WORK/r.log" 2>&1
  if diff -q "$c.clstr" "$r.clstr" >/dev/null 2>&1 && diff -q "$c" "$r" >/dev/null 2>&1; then
    pass=$((pass+1))
  else
    echo "FAIL $label"
    diff "$c.clstr" "$r.clstr" 2>&1 | head -8
    diff -q "$c" "$r" 2>&1 | head -2
    fail=$((fail+1)); failed+=("$label")
  fi
}
est2d() { run_2d cd-hit-est-2d "$@"; }
aa2d()  { run_2d cd-hit-2d     "$@"; }

# cd-hit-div splits a database into `-div` files
run_div() {                       # run_div <input> <label> <div>
  local input="$1" label="$2" div="$3"
  local c="$WORK/c.$label" r="$WORK/r.$label"
  "$CPP/cd-hit-div" -i "$input" -o "$c" -div "$div" > "$WORK/c.log" 2>&1
  "$RS/cd-hit-div"  -i "$input" -o "$r" -div "$div" > "$WORK/r.log" 2>&1
  local ok=1
  for k in $(seq 0 $((div-1))); do
    [ -f "$c-$k" ] || continue
    diff -q "$c-$k" "$r-$k" >/dev/null 2>&1 || ok=0
  done
  # the Rust side must not produce extra or missing parts either
  local nc nr
  nc=$(ls "$c"-* 2>/dev/null | wc -l); nr=$(ls "$r"-* 2>/dev/null | wc -l)
  [ "$nc" = "$nr" ] || ok=0
  if [ $ok = 1 ] && [ "$nc" != "0" ]; then
    pass=$((pass+1))
  else
    echo "FAIL $label (parts: $nc vs $nr)"
    fail=$((fail+1)); failed+=("$label")
  fi
}

# ---- cd-hit-est -----------------------------------------------------------
for c in 0.80 0.90 0.95 0.99 1.00; do
  est "$DATA/dna_small.fna" "est_small_c$c"  -c $c -T 1
  est "$DATA/dna_mixed.fna" "est_mixed_c$c"  -c $c -T 1
done
for n in 4 5 6 8 9 10 11; do
  est "$DATA/dna_mixed.fna" "est_mixed_n$n" -c 0.9 -n $n -T 1
done
est "$DATA/dna_mixed.fna" "est_norev"     -c 0.9 -r 0 -T 1
est "$DATA/dna_mixed.fna" "est_best"      -c 0.9 -g 1 -T 1
est "$DATA/dna_mixed.fna" "est_print"     -c 0.9 -p 1 -T 1
est "$DATA/dna_mixed.fna" "est_desc0"     -c 0.9 -d 0 -T 1
est "$DATA/dna_mixed.fna" "est_desc100"   -c 0.9 -d 100 -T 1
est "$DATA/dna_mixed.fna" "est_local"     -c 0.9 -G 0 -aS 0.8 -T 1
est "$DATA/dna_mixed.fna" "est_cov"       -c 0.9 -aL 0.7 -aS 0.7 -T 1
est "$DATA/dna_mixed.fna" "est_covctl"    -c 0.9 -AL 40 -AS 40 -T 1
est "$DATA/dna_mixed.fna" "est_lendiff"   -c 0.9 -s 0.8 -S 100 -T 1
est "$DATA/dna_mixed.fna" "est_band"      -c 0.9 -b 5 -T 1
est "$DATA/dna_mixed.fna" "est_band100"   -c 0.9 -b 100 -T 1
est "$DATA/dna_mixed.fna" "est_minlen"    -c 0.9 -l 200 -T 1
est "$DATA/dna_mixed.fna" "est_mask"      -c 0.9 -mask N -T 1
est "$DATA/dna_mixed.fna" "est_scores"    -c 0.9 -match 3 -mismatch -3 -T 1
est "$DATA/dna_mixed.fna" "est_gap"       -c 0.9 -gap -4 -gap-ext -1 -T 1
est "$DATA/dna_mixed.fna" "est_sortc"     -c 0.9 -sc 1 -T 1
est "$DATA/dna_mixed.fna" "est_unmatch"   -c 0.9 -uL 0.3 -uS 0.3 -U 50 -T 1
est "$DATA/dna_mixed.fna" "est_mem200"    -c 0.9 -M 200 -T 1
est "$DATA/dna_mixed.fna" "est_mem0"      -c 0.9 -M 0 -T 1
est "$DATA/dna_mixed.fna.gz" "est_gz"     -c 0.9 -T 1
est "$DATA/dna_dup.fna"   "est_dup"       -c 0.98 -T 1
est "$DATA/dna_split.fna" "est_split"     -c 0.90 -T 1
est "$DATA/dna_dup.fna"   "est_dup_trim"  -c 0.95 -cx 100 -T 1
est "$DATA/dna_small.fna" "est_dist"      -D 0.05 -T 1

# ---- long N50: ComputeTableLimits overflows its `int` expression ----------
# `500*max_len + 500000*typical_len + 50000000` is int arithmetic upstream and wraps
# negative once the N50 passes ~4294. `(size_t)` of that negative double is ~1.8e19, i.e.
# effectively no limit; a saturating conversion would give 0 and spin the block loop in
# DoClustering forever. Check with: the run prints
# "Max number of word counting entries: 18446744072467088320".
for c in 0.90 0.95 0.99; do
  est "$DATA/dna_longn50.fna" "est_longn50_c$c" -c $c -M 0 -T 1
done
est "$DATA/dna_longn50.fna" "est_longn50_mem"  -c 0.95 -M 200 -T 1
est "$DATA/dna_longn50.fna" "est_longn50_n8"   -c 0.95 -n 8 -M 0 -T 1
est "$DATA/dna_longn50.fna" "est_longn50_dist" -D 0.05 -M 0 -T 1
p454 "$DATA/dna_longn50.fna" "454_longn50"     -c 0.98 -M 0 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_longn50.fna" "est2d_longn50" -c 0.9 -M 0 -T 1
# Guard the guard: if the N50 ever drops back under the overflow point these cases stop
# testing anything.
entries=$("$RS/cd-hit-est" -i "$DATA/dna_longn50.fna" -o "$WORK/n50" -c 0.95 -M 0 -T 1 2>&1 |
            sed -n 's/^Max number of word counting entries: //p')
# (a bare `-gt` would overflow bash's own signed 64-bit arithmetic here)
if [ "${#entries}" -ge 19 ]; then
  pass=$((pass+1))
else
  echo "FAIL longn50_guard: max_entries=$entries, expected the overflowed ~1.8e19"
  fail=$((fail+1)); failed+=("longn50_guard")
fi

# ---- cd-hit-2d / cd-hit-est-2d -------------------------------------------
for c in 0.80 0.90 0.95 0.99; do
  est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_c$c" -c $c -T 1
done
est2d "$DATA/dna_mixed.fna" "$DATA/dna_small.fna" "est2d_swap"   -c 0.9 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_s2"     -c 0.9 -s2 0.8 -S2 50 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_norev"  -c 0.9 -r 0 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_best"   -c 0.9 -g 1 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_print"  -c 0.9 -p 1 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_cov"    -c 0.9 -aL 0.7 -aS 0.7 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_mem"    -c 0.9 -M 100 -T 1
est2d "$DATA/dna_small.fna" "$DATA/dna_mixed.fna" "est2d_n8"     -c 0.9 -n 8 -T 1
est2d "$DATA/dna_dup.fna"   "$DATA/dna_small.fna" "est2d_dup"    -c 0.95 -T 1

for c in 0.70 0.90 0.98; do
  aa2d "$DATA/prot_small.faa" "$DATA/prot_mixed.faa" "aa2d_c$c" -c $c -T 1
done
aa2d "$DATA/prot_mixed.faa" "$DATA/prot_small.faa" "aa2d_swap"  -c 0.9 -T 1
aa2d "$DATA/prot_small.faa" "$DATA/prot_mixed.faa" "aa2d_s2"    -c 0.9 -s2 0.8 -S2 40 -T 1
aa2d "$DATA/prot_small.faa" "$DATA/prot_mixed.faa" "aa2d_best"  -c 0.9 -g 1 -T 1
aa2d "$DATA/prot_small.faa" "$DATA/prot_mixed.faa" "aa2d_print" -c 0.9 -p 1 -T 1
aa2d "$DATA/prot_small.faa" "$DATA/prot_mixed.faa" "aa2d_cov"   -c 0.9 -aL 0.7 -aS 0.7 -T 1
aa2d "$DATA/prot_small.faa" "$DATA/prot_mixed.faa" "aa2d_n4"    -c 0.9 -n 4 -T 1
aa2d "$DATA/prot_small.faa" "$DATA/prot_mixed.faa" "aa2d_t0"    -c 0.9 -t 0 -T 1

# ---- cd-hit-454 -----------------------------------------------------------
p454 "$DATA/dna_dup.fna"   "454_dup"    -c 0.98 -T 1
p454 "$DATA/dna_small.fna" "454_small"  -c 0.98 -T 1
p454 "$DATA/dna_mixed.fna" "454_mixed"  -c 0.98 -T 1
p454 "$DATA/dna_dup.fna"   "454_indel"  -c 0.98 -D 3 -T 1
p454 "$DATA/dna_dup.fna"   "454_c99"    -c 0.99 -T 1

# ---- cd-hit-div -----------------------------------------------------------
run_div "$DATA/dna_mixed.fna"  "div_dna_4"  4
run_div "$DATA/prot_mixed.faa" "div_aa_3"   3
run_div "$DATA/dna_small.fna"  "div_dna_10" 10

# ---- cd-hit ---------------------------------------------------------------
for c in 0.50 0.70 0.90 0.98 1.00; do
  prot "$DATA/prot_small.faa" "aa_small_c$c" -c $c -T 1
  prot "$DATA/prot_mixed.faa" "aa_mixed_c$c" -c $c -T 1
done
for n in 2 3 4 5; do
  prot "$DATA/prot_mixed.faa" "aa_mixed_n$n" -c 0.9 -n $n -T 1
done
for t in 0 1 2 3 4 5; do
  prot "$DATA/prot_mixed.faa" "aa_mixed_t$t" -c 0.9 -t $t -T 1
done
prot "$DATA/prot_mixed.faa" "aa_best"     -c 0.9 -g 1 -T 1
prot "$DATA/prot_mixed.faa" "aa_print"    -c 0.9 -p 1 -T 1
prot "$DATA/prot_mixed.faa" "aa_local"    -c 0.9 -G 0 -aS 0.8 -T 1
prot "$DATA/prot_mixed.faa" "aa_cov"      -c 0.9 -aL 0.7 -aS 0.7 -T 1
prot "$DATA/prot_mixed.faa" "aa_lendiff"  -c 0.9 -s 0.8 -S 60 -T 1
prot "$DATA/prot_mixed.faa" "aa_band"     -c 0.9 -b 5 -T 1
prot "$DATA/prot_mixed.faa" "aa_desc0"    -c 0.9 -d 0 -T 1
prot "$DATA/prot_mixed.faa" "aa_sortc"    -c 0.9 -sc 1 -T 1
prot "$DATA/prot_mixed.faa" "aa_mem200"   -c 0.9 -M 200 -T 1
prot "$DATA/prot_mixed.faa" "aa_mem0"     -c 0.9 -M 0 -T 1
prot "$DATA/prot_split.faa" "aa_split"   -c 0.90 -T 1
prot "$DATA/prot_mixed.faa" "aa_unmatch"  -c 0.9 -uL 0.3 -uS 0.3 -T 1
prot "$DATA/prot_small.faa" "aa_dist"     -D 0.1 -T 1

# ---- coverage / length options driven outside [0,1] -----------------------
# Panaroo runs cd-hit with -aS/-AS/-AL/-S far out of range. `ControlShortCoverage` then
# evaluates `(int)(len * 99999999.0)`, which does not fit in an int: on x86-64 the C++
# lands on INT_MIN and the very next line overwrites it, while a saturating conversion
# would keep INT_MAX and reject every alignment. Same story for -aL and for a tiny -s in
# `upper_bound_length_rep`.
PANAROO="-s None -aL 0.0 -AL 99999999 -aS 99999999 -AS 99999999 -M 0 -d 999 -g 1 -n 2"
for c in 0.70 0.80 0.90 0.95 0.99; do
  prot "$DATA/prot_pair.faa"    "aa_panaroo_pair_c$c" -c $c $PANAROO -T 1
  prot "$DATA/prot_panaroo.faa" "aa_panaroo_c$c"      -c $c $PANAROO -T 1
done
prot "$DATA/prot_mixed.faa" "aa_aS_huge"  -c 0.9 -aS 99999999 -AS 99999999 -T 1
prot "$DATA/prot_mixed.faa" "aa_aL_huge"  -c 0.9 -aL 99999999 -AL 99999999 -T 1
prot "$DATA/prot_mixed.faa" "aa_G0_aS_huge" -c 0.9 -G 0 -aS 99999999 -AS 99999999 -T 1
prot "$DATA/prot_mixed.faa" "aa_s_tiny"   -c 0.9 -s 0.000000001 -T 1
prot "$DATA/prot_mixed.faa" "aa_aL_tiny"  -c 0.9 -aL 0.000000001 -T 1
# `upper_bound_length_rep` divides in double precision, not float: `-s 0.001` on a
# 333aa sequence must give 333000, not 332999.
for sm in 0.001 0.002 0.008 0.016; do
  prot "$DATA/prot_mixed.faa" "aa_s_$sm"  -c 0.9 -s $sm -T 1
  est  "$DATA/dna_mixed.fna"  "est_s_$sm" -c 0.9 -s $sm -T 1
done
est  "$DATA/dna_mixed.fna"  "est_aS_huge" -c 0.9 -aS 99999999 -AS 99999999 -T 1
est  "$DATA/dna_mixed.fna"  "est_aL_huge" -c 0.9 -aL 99999999 -AL 99999999 -T 1
est  "$DATA/dna_mixed.fna"  "est_s_tiny"  -c 0.9 -s 0.000000001 -T 1

# ---- -D distance mode: `ComputeRequiredBases` takes its own branch ---------
for d in 0.01 0.02 0.05 0.10 0.20; do
  prot "$DATA/prot_small.faa" "aa_dist_$d"   -D $d -T 1
  prot "$DATA/prot_mixed.faa" "aa_distm_$d"  -D $d -T 1
  est  "$DATA/dna_small.fna"  "est_dist_$d"  -D $d -T 1
  est  "$DATA/edge_short.fna" "est_dists_$d" -D $d -T 1
done
prot "$DATA/prot_small.faa" "aa_dist_n2"  -D 0.05 -n 2 -T 1
prot "$DATA/prot_small.faa" "aa_dist_g1"  -D 0.05 -g 1 -T 1

# ---- awkward inputs -------------------------------------------------------
for f in empty one no_trailing_nl crlf blank_lines lowercase dup_names long_header \
         no_space_header tab_header space_after_gt short one_per_line iupac bad_chars \
         empty_seq identical homopolymer all_n; do
  est "$DATA/edge_$f.fna" "edge_est_$f" -c 0.9 -T 1
  est "$DATA/edge_$f.fna" "edge_est_${f}_d0" -c 0.9 -d 0 -T 1
done
prot "$DATA/edge_stops.faa" "edge_aa_stops" -c 0.9 -T 1
prot "$DATA/edge_bad_chars.fna" "edge_aa_bad" -c 0.9 -T 1

# ---- refusals: both sides must reject the same way ------------------------
prot "$DATA/prot_mixed.faa" "aa_thd_too_low"  -c 0.5 -T 1
prot "$DATA/prot_mixed.faa" "aa_bad_word"     -c 0.9 -n 9 -T 1
est  "$DATA/dna_mixed.fna"  "est_thd_too_low" -c 0.5 -T 1
est  "$DATA/dna_mixed.fna"  "est_bad_band"    -c 0.9 -b 0 -T 1

# ---- paired end -----------------------------------------------------------
pe_case() {
  local label="$1"; shift
  local c="$WORK/c.$label" r="$WORK/r.$label"
  "$CPP/cd-hit-est" -i "$DATA/pe_R1.fq" -j "$DATA/pe_R2.fq" -o "$c" -op "$c.pe" -P 1 "$@" \
      > "$WORK/c.log" 2>&1
  "$RS/cd-hit-est"  -i "$DATA/pe_R1.fq" -j "$DATA/pe_R2.fq" -o "$r" -op "$r.pe" -P 1 "$@" \
      > "$WORK/r.log" 2>&1
  if diff -q "$c.clstr" "$r.clstr" >/dev/null 2>&1 \
     && diff -q "$c" "$r" >/dev/null 2>&1 \
     && diff -q "$c.pe" "$r.pe" >/dev/null 2>&1; then
    pass=$((pass+1))
  else
    echo "FAIL $label"
    diff "$c.clstr" "$r.clstr" 2>&1 | head -8
    fail=$((fail+1)); failed+=("$label")
  fi
}
pe_case pe_basic  -c 0.95 -T 1
pe_case pe_trim   -c 0.95 -cx 80 -cy 80 -T 1
pe_case pe_sortf  -c 0.95 -sf 1 -T 1

# ---- threading: -T must not change the answer -----------------------------
# The C++ OpenMP path is racy and gives a different clustering at each -T (and between
# runs at the same -T), so the reference for every thread count is the C++ -T 1 output.
thread_case() {                   # thread_case <prog> <input> <label> <opts...>
  local prog="$1" input="$2" label="$3"; shift 3
  local c="$WORK/c.$label"
  "$CPP/$prog" -i "$input" -o "$c" -T 1 "$@" > "$WORK/c.log" 2>&1
  for t in 2 4 8; do
    local r="$WORK/r.$label.$t"
    "$RS/$prog" -i "$input" -o "$r" -T $t "$@" > "$WORK/r.log" 2>&1
    if diff -q "$c.clstr" "$r.clstr" >/dev/null 2>&1 && diff -q "$c" "$r" >/dev/null 2>&1; then
      pass=$((pass+1))
    else
      echo "FAIL $label -T $t (differs from C++ -T 1)"
      diff "$c.clstr" "$r.clstr" 2>&1 | head -6
      fail=$((fail+1)); failed+=("$label-T$t")
    fi
  done
}
thread_case cd-hit-est "$DATA/dna_mixed.fna" "thr_est"      -c 0.90
thread_case cd-hit-est "$DATA/dna_dup.fna"   "thr_est_dup"  -c 0.95
thread_case cd-hit     "$DATA/prot_mixed.faa" "thr_aa"      -c 0.90
thread_case cd-hit-454 "$DATA/dna_dup.fna"   "thr_454"      -c 0.98
thread_case cd-hit-est "$DATA/dna_longn50.fna" "thr_longn50" -c 0.95 -M 0

# The cases above are all under 1000 sequences, and `DoClustering`'s -T block split only
# engages above that (`if( m > i + 1E3 ) ...`). They therefore still run as a single block
# at every -T, and the parallel screening phase gets an empty range. The two inputs below
# are over 2000 sequences, so -T really does cut the block and the workers really do run.
#
# Under -g 1 the split is invisible: `CheckOne` never marks a sequence redundant, it only
# records the best hit so far in `seq.identity` and skips any candidate that cannot beat
# it, so the winner is the global best whichever block it was met in. Hence the strict
# comparison against the C++ -T 1 output still applies.
thread_case cd-hit     "$DATA/prot_split.faa" "thr_aa_split_g"   -c 0.90 -g 1
thread_case cd-hit     "$DATA/prot_split.faa" "thr_aa_split99_g" -c 0.99 -n 2 -g 1
thread_case cd-hit-est "$DATA/dna_split.fna"  "thr_est_split_g"  -c 0.90 -g 1
thread_case cd-hit-est "$DATA/dna_split.fna"  "thr_est_split99_g" -c 0.99 -n 8 -g 1

# Under -g 0 it is not: `CheckOne` stops at its *first* hit, and the candidate order
# depends on which representatives share a word with the query, so meeting an earlier
# block's representatives first can pick a different -- equally valid -- representative.
# Upstream's -T has exactly the same sensitivity (`cdhit/cd-hit -T 2` on prot_split.faa
# differs from its own -T 1), so there is no fixed answer to compare against. What can be
# required, and is more than upstream manages, is that each -T is reproducible: the same
# input and the same -T must give the same clustering every time.
repeat_case() {                   # repeat_case <prog> <input> <label> <opts...>
  local prog="$1" input="$2" label="$3"; shift 3
  for t in 2 4 8; do
    local a="$WORK/rep.$label.$t.a" b="$WORK/rep.$label.$t.b"
    "$RS/$prog" -i "$input" -o "$a" -T $t "$@" > "$WORK/r.log" 2>&1
    "$RS/$prog" -i "$input" -o "$b" -T $t "$@" > "$WORK/r.log" 2>&1
    if diff -q "$a.clstr" "$b.clstr" >/dev/null 2>&1 && diff -q "$a" "$b" >/dev/null 2>&1; then
      pass=$((pass+1))
    else
      echo "FAIL $label -T $t (two runs at the same -T disagree)"
      diff "$a.clstr" "$b.clstr" 2>&1 | head -6
      fail=$((fail+1)); failed+=("$label-T$t")
    fi
  done
}
repeat_case cd-hit     "$DATA/prot_split.faa" "rep_aa_split"  -c 0.90
repeat_case cd-hit-est "$DATA/dna_split.fna"  "rep_est_split" -c 0.90

# The -M route to several blocks, which is what a memory-limited run takes and which the
# -T split has to compose with. Squeeze -M just above the memory floor.
thread_case cd-hit-est "$DATA/dna_mixed.fna"  "thr_est_mb"  -c 0.90 -n 8 -M 14
thread_case cd-hit-est "$DATA/dna_dup.fna"    "thr_dup_mb"  -c 0.95 -n 8 -M 14
thread_case cd-hit     "$DATA/prot_mixed.faa" "thr_aa_mb"   -c 0.90 -n 4 -M 14
thread_case cd-hit-est "$DATA/dna_small.fna"  "thr_small_mb" -c 0.95 -n 8 -M 14

# Guard the guard: if a change to the memory accounting moves the block boundaries, the
# -M cases would silently go back to being single-block and stop testing anything. The
# -T cases are guarded the same way, at -T 8 -- a regression that made -T a no-op again
# would leave them at one block, and every diff above would still pass.
for spec in "cd-hit-est|$DATA/dna_mixed.fna|-c 0.90 -n 8 -M 14|1" \
            "cd-hit-est|$DATA/dna_dup.fna|-c 0.95 -n 8 -M 14|1" \
            "cd-hit|$DATA/prot_mixed.faa|-c 0.90 -n 4 -M 14|1" \
            "cd-hit-est|$DATA/dna_small.fna|-c 0.95 -n 8 -M 14|1" \
            "cd-hit|$DATA/prot_split.faa|-c 0.90|8" \
            "cd-hit-est|$DATA/dna_split.fna|-c 0.90|8"; do
  prog="${spec%%|*}"; rest="${spec#*|}"; input="${rest%%|*}"; rest="${rest#*|}"
  opts="${rest%%|*}"; t="${rest##*|}"
  blocks=$("$RS/$prog" -i "$input" -o "$WORK/blk" $opts -T $t 2>&1 |
             grep -c "comparing sequences from")
  if [ "$blocks" -gt 1 ]; then
    pass=$((pass+1))
  else
    echo "FAIL multiblock_guard $prog $(basename "$input") -T $t: $blocks block(s), expected >1"
    echo "     the -T cases for this input are no longer exercising the parallel path"
    fail=$((fail+1)); failed+=("multiblock_guard_${prog}_T$t")
  fi
done

# 2D under threads, likewise against the C++ -T 1 output
for t in 2 4 8; do
  c="$WORK/c.thr2d"; r="$WORK/r.thr2d.$t"
  "$CPP/cd-hit-est-2d" -i "$DATA/dna_small.fna" -i2 "$DATA/dna_mixed.fna" -o "$c" -c 0.9 -T 1       > "$WORK/c.log" 2>&1
  "$RS/cd-hit-est-2d"  -i "$DATA/dna_small.fna" -i2 "$DATA/dna_mixed.fna" -o "$r" -c 0.9 -T $t       > "$WORK/r.log" 2>&1
  if diff -q "$c.clstr" "$r.clstr" >/dev/null 2>&1 && diff -q "$c" "$r" >/dev/null 2>&1; then
    pass=$((pass+1))
  else
    echo "FAIL thr2d -T $t"; fail=$((fail+1)); failed+=("thr2d-T$t")
  fi
done

echo
echo "passed: $pass   failed: $fail"
if [ $fail -gt 0 ]; then printf 'failing: %s\n' "${failed[*]}"; exit 1; fi
