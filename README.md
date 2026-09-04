# cdhit-rs

A Rust translation of [CD-HIT](https://github.com/weizhongli/cdhit) 4.8.1. 

Note that the original code is nondeterministic when using multiple threads: its `.clstr`
differs between two runs at the same `-T`, and differs again from `-T 1`. This port is
deterministic — two runs at the same `-T` always agree — and under `-g 1` its output is
bit-identical at every thread count and matches the C++ `-T 1`. See [Threading](#threading).

* 2026-09-03: First translation. There are some intricate details in the code so more testing needed


## CLI

All six C++ programs are translated and available as subcommands:

```
cdhit est -i in.fna -o out -c 0.95        # the short form
cdhit cd-hit-est -i in.fna -o out -c 0.95 # the upstream name also works
ln -s cdhit cd-hit-est && ./cd-hit-est …  # and so does a symlink, for existing scripts
```

| command | short form | from | what it does |
| --- | --- | --- | --- |
| `cd-hit` | `prot` | `cdhit.c++` | protein clustering |
| `cd-hit-est` | `est` | `cdhit-est.c++` | nucleotide clustering, both strands |
| `cd-hit-2d` | `2d` | `cdhit-2d.c++` | compare db2 against db1, db1 unclustered |
| `cd-hit-est-2d` | `est-2d` | `cdhit-est-2d.c++` | the same for nucleotides |
| `cd-hit-454` | `454` | `cdhit-454.c++` | 454 duplicate removal, indel-tolerant |
| `cd-hit-div` | `div` | `cdhit-div.c++` | split a database into N parts |

Out of scope: `cd-hit-auxtools` (a separate codebase), `psi-cd-hit` (needs BLAST), and the
Perl post-processing scripts.

## Status

The translation is faithful. Across every input and option combination in `tests/` — 263
cases — the `.clstr` file and the representative FASTA are **byte-identical** to the C++
binary's, and a randomised differential harness adds several hundred more.

It is also faster everywhere measured — 1.1x to 2.4x depending on configuration. Best of
three, user CPU time, `-T 1`, on an otherwise-quiet Xeon Gold 6138; every case produced
byte-identical output to the C++:

| case | C++ | Rust | speedup (C++ / Rust) |
| --- | --- | --- | --- |
| `cd-hit-est -n 5` | 15.74s | 14.02s | 1.12x |
| `cd-hit-est -n 8` | 2.76s | 2.22s | 1.24x |
| `cd-hit-est -n 10` | 2.47s | 1.37s | 1.80x |
| `cd-hit-est -n 11` | 2.32s | 1.46s | 1.59x |
| `cd-hit-est -c 0.90` | 94.26s | 50.44s | 1.87x |
| `cd-hit -n 5` | 4.13s | 1.75s | 2.36x |
| `cd-hit -n 4` | 5.91s | 4.26s | 1.39x |
| `cd-hit -c 0.70` | 4.35s | 2.63s | 1.65x |
| `cd-hit-454` | 1.39s | 0.57s | 2.44x |
| `cd-hit-est-2d` | 5.34s | 2.91s | 1.84x |
| `cd-hit-2d` | 9.68s | 5.35s | 1.81x |

Reproduce with `BENCH_DNA=<dna.fna> BENCH_AA=<prot.faa> ./tests/bench.sh`.

Peak RSS is slightly lower than the C++ too (89 MB vs 93 MB on the `-n 10` case), and grows
far less with `-T`: 21 MB at `-T 1` and 11 MB at `-T 20` on the protein case in
[Threading](#threading), against 32 MB and 344 MB for the C++.


```
cargo build --release
cargo test                   # unit tests for the leaf functions
./tests/difftest.sh          # differential tests against the C++ binaries
./tests/fuzz.py              # randomised differential test, aimed at the SIMD kernel
./tests/bench.sh             # timing comparison
```

`tests/fuzz.py [rust-binary] [cases]` generates random nucleotide inputs — varying lengths,
N and masked-base placement, homopolymer runs, identity threshold, word length and band
width — and requires byte-identical output from both binaries on each. `SEED=n` picks the
stream. It exists because the vectorised diagonal kernel is `unsafe` and its correctness
depends on exactly those properties; a fixed corpus does not vary them enough.

`difftest.sh` covers all six programs: identity thresholds, word lengths, coverage and
length-difference options, `-g`/`-G`/`-p`/`-d`/`-sc`/`-sf`, custom scoring, masking,
distance mode, gzip input, FASTQ, paired-end, memory limits, thread counts, and a set of
awkward inputs (empty files, CRLF, missing trailing newline, lowercase, all-N, homopolymers,
IUPAC codes, invalid characters, headers with no whitespace). Point it at the binaries with
`CPP_DIR` / `RS_DIR`; regenerate the inputs with `python3 tests/gen_data.py && python3
tests/gen_edge.py`.


## Threading

`SequenceDB::DoClustering` splits the database into blocks: a block is clustered serially
into a word table, then every later sequence is screened against that finished table, which
is the part that parallelises. Upstream cuts the block down to `(N - i) / (2 + T)` once
`-T` is in play; without that cut, a database that fits one word table is a single block
and the screening phase gets an empty range. This port applies the cut, from upstream's own
formula (`cdhit-common.c++:3098`).

Best of five, 5,359 protein representatives, `-c 0.99 -g 1 -n 2 -M 0`, 40-core Xeon Gold
6138. All twenty runs produced byte-identical output:

| `-T` | C++ wall / CPU | Rust wall / CPU |
| --- | --- | --- |
| 1 | 2.96s / 2.93s | 2.49s / 2.47s |
| 4 | 1.13s / 3.61s | 1.05s / 2.90s |
| 8 | 0.70s / 3.99s | 0.81s / 3.81s |
| 20 | 0.58s / 6.04s | 0.44s / 4.82s |

On a redundancy-heavy input the split also *removes* work: a sequence matched by an early
block's table is skipped by every later block, where a single block would have compared it
against everything. On 21,436 protein sequences at `-c 0.98 -g 1 -n 2` (four mutated copies
of the input above), CPU time falls from 6.7s at `-T 1` to 1.9s at `-T 20` and wall time to
1.3s; the C++ manages 2.5s at `-T 8` and gets slower again at `-T 20`.

Determinism: two runs at the same `-T` always agree, unlike upstream, whose `may_stop` /
`self_stop` / `stop` flags are written and read by several threads under nothing stronger
than `#pragma omp flush`. Across thread counts:

* **`-g 1`**: `-T` cannot change the answer. `CheckOne` never marks a sequence redundant in
  that mode — it records the best hit so far in `seq.identity` and skips any candidate that
  cannot beat it — so the winner is the global best whichever block it was met in.
* **`-g 0`** (the default): `CheckOne` stops at its *first* hit, so meeting an earlier
  block's representatives first can select a different, equally valid representative.
  Upstream's `-T` has the same sensitivity. It rarely bites — on both the protein input
  above and a 5,359-sequence nucleotide one every `-T` still gives byte-identical output to
  `-T 1` and to the C++ — but `tests/data/prot_split.faa` is deliberately noisy enough to
  show it, for the C++ as well as here.

Upstream's parallel driver also overlaps the serial block build with the screening, out of
a second word table. That is not reproduced, and is the remaining headroom: the serial block
loop is about 6% of the comparisons at `-T 20` and 9% at `-T 8`, which is where the C++
`-T 8` figure above comes from. See `NAMING.md`.

## Licence

GNU General Public License, version 2

## Citation

Please cite the original software:

1. Weizhong Li & Adam Godzik. CD-HIT: a fast program for clustering and comparing large sets
   of protein or nucleotide sequences. Bioinformatics (2006) 22:1658-1659.
2. Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. CD-HIT: accelerated for
   clustering the next generation sequencing data. Bioinformatics (2012) 28:3150-3152.
3. Beifang Niu, Limin Fu, Shulei Sun & Weizhong Li. Artificial and natural duplicates in
   pyrosequencing reads of metagenomic data. BMC Bioinformatics (2010) 11:187.

If you use our translation, we recommend that you also cite the precise version you use. If you link to [crates.io](http://crates.io), you can cite the version number;
but if you link to our Git repository, for reproducibility, it is better that you provide the URL to the repository and the git hash (Github lists it high up on the page as 7 letters, under the Code button, e.g. '21751cd')

In addition, we appreciate if you cite the paper below describing the translation approach. If for some reason you struggle with journal citation limits, please prioritizing citing the original software over our translation paper.

> Johan Henriksson. Static analysis-guided agentic AI translation enables Rust as a full stack bioinformatics language. arXiv:2608.13029, 2026. https://doi.org/10.48550/arXiv.2608.13029

