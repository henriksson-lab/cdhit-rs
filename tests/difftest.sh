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
est "$DATA/dna_dup.fna"   "est_dup_trim"  -c 0.95 -cx 100 -T 1
est "$DATA/dna_small.fna" "est_dist"      -D 0.05 -T 1

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
prot "$DATA/prot_mixed.faa" "aa_unmatch"  -c 0.9 -uL 0.3 -uS 0.3 -T 1
prot "$DATA/prot_small.faa" "aa_dist"     -D 0.1 -T 1

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
