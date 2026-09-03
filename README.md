# cdhit-rs

A Rust translation of [CD-HIT](https://github.com/weizhongli/cdhit) 4.8.1. The C++ original
is vendored at `cdhit/` and is the reference this port is checked against.

All six C++ programs are translated:

| binary | from | what it does |
| --- | --- | --- |
| `cd-hit` | `cdhit.c++` | protein clustering |
| `cd-hit-est` | `cdhit-est.c++` | nucleotide clustering, both strands |
| `cd-hit-2d` | `cdhit-2d.c++` | compare db2 against db1, db1 unclustered |
| `cd-hit-est-2d` | `cdhit-est-2d.c++` | the same for nucleotides |
| `cd-hit-454` | `cdhit-454.c++` | 454 duplicate removal, indel-tolerant |
| `cd-hit-div` | `cdhit-div.c++` | split a database into N parts |

Out of scope: `cd-hit-auxtools` (a separate codebase), `psi-cd-hit` (needs BLAST), and the
Perl post-processing scripts.

## Status

The translation is faithful. Across every input and option combination in `tests/` — 165
cases — the `.clstr` file and the representative FASTA are **byte-identical** to the C++
binary's.

It is also faster everywhere measured. User CPU, best of two, on a shared Xeon Gold 6138;
`same` means the output matched byte for byte:

| case | C++ | Rust | |
| --- | --- | --- | --- |
| `cd-hit-est -n 5` | 15.77s | 13.76s | 1.15x |
| `cd-hit-est -n 8` | 2.77s | 2.51s | 1.10x |
| `cd-hit-est -n 10` | 2.43s | 1.81s | 1.34x |
| `cd-hit-est -n 11` | 2.32s | 1.61s | 1.44x |
| `cd-hit-est -c 0.90` | 92.90s | 86.67s | 1.07x |
| `cd-hit -n 5` | 4.07s | 2.77s | 1.47x |
| `cd-hit -n 4` | 5.83s | 4.83s | 1.21x |
| `cd-hit -c 0.70` | 4.33s | 3.71s | 1.17x |
| `cd-hit-454` | 1.40s | 0.74s | 1.89x |

Peak RSS is slightly lower than the C++ too (89 MB vs 93 MB on the `-n 10` case).

```
cargo build --release
cargo test                   # unit tests for the leaf functions
./tests/difftest.sh          # differential tests against the C++ binaries
./tests/bench.sh             # timing comparison
```

`difftest.sh` covers all six programs: identity thresholds, word lengths, coverage and
length-difference options, `-g`/`-G`/`-p`/`-d`/`-sc`/`-sf`, custom scoring, masking,
distance mode, gzip input, FASTQ, paired-end, memory limits, thread counts, and a set of
awkward inputs (empty files, CRLF, missing trailing newline, lowercase, all-N, homopolymers,
IUPAC codes, invalid characters, headers with no whitespace). Point it at the binaries with
`CPP_DIR` / `RS_DIR`; regenerate the inputs with `python3 tests/gen_data.py && python3
tests/gen_edge.py`.

## Layout

```
Cargo.toml          one crate, `cdhit-rs`; the six front-ends come from src/bin/
build.rs            stamps the build date into the banner, as __DATE__ does upstream
src/                the translation
tests/              unit tests, the differential suite, and its generated inputs
cdhit/              the vendored C++ 4.8.1 this is checked against
```

| file | translated from |
| --- | --- |
| `src/cdhit_common.rs` | `cdhit-common.h` + `cdhit-common.c++` |
| `src/cdhit_utility.rs` | `cdhit-utility.h` + `cdhit-utility.c++` |
| `src/bin/*.rs` | the corresponding `cdhit-*.c++` main |

One Rust module per C++ file pair, items in their original order, one Rust function per C++
function. `NAMING.md` has the C++ -> Rust name mapping, the deliberate structural changes,
and the upstream defects that are reproduced on purpose.

## Threading

`-T` is a pure speed knob here: the output is bit-identical at every thread count, and
matches the C++ `-T 1` output. The upstream OpenMP path is racy — it produces a different
clustering from run to run at the same `-T` — so there was no fixed behaviour to be
compatible with. See `NAMING.md`.
