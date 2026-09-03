# cdhit-rs

A Rust translation of [CD-HIT](https://github.com/weizhongli/cdhit) 4.8.1. 

Note that the original code is nondeterministic when using multiple threads: its `.clstr`
differs between two runs at the same `-T`, and differs again from `-T 1`. This port is
deterministic — `-T` changes only speed, and the output matches the C++ `-T 1` at every
thread count.

* 2026-09-03: First translation. There are some intricate details in the code so more testing needed


## CLI

All six C++ programs are translated. Upstream ships them as six executables; here they are
subcommands of one binary, which is the only change to how you invoke them:

```
cdhit est -i in.fna -o out -c 0.95        # the short form
cdhit cd-hit-est -i in.fna -o out -c 0.95 # the upstream name also works
ln -s cdhit cd-hit-est && ./cd-hit-est …  # and so does a symlink, for existing scripts
```

`cdhit` on its own lists the commands; `cdhit <command> -h` gives that command's options.
Each command is the faithful translation of the corresponding `main`, unchanged.


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

The translation is faithful. Across every input and option combination in `tests/` — 241
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

Peak RSS is slightly lower than the C++ too (89 MB vs 93 MB on the `-n 10` case), and no
longer grows with `-T`.


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

`-T` is a pure speed knob here: the output is bit-identical at every thread count, and
matches the C++ `-T 1` output. The upstream OpenMP path is racy — it produces a different
clustering from run to run at the same `-T` — so there was no fixed behaviour to be
compatible with. See `NAMING.md`.

## Licence

CD-HIT 4.8.1 is under the GNU General Public License, version 2. `cdhit/license.txt` is the
GPLv2 text verbatim, and the user's guide (`cdhit/doc/cdhit-user-guide.wiki`, "Installation")
says "It can be copied under the GNU General Public License version 2 (GPLv2)". Neither adds
an "or any later version" clause.

GPLv2 section 0 counts a translation as a modification — "either verbatim or with
modifications and/or translated into another language" — so this port is a derivative work of
CD-HIT and is distributed under the same licence, GPLv2. It cannot be relicensed. `LICENSE`
at the root is a copy of `cdhit/license.txt`.

CD-HIT was written by Weizhong Li (UCSD, San Diego Supercomputer Center) in Adam Godzik's lab
at the Burnham Institute, and modified by Limin Fu (Center for Research in Biological Systems,
UCSD); the file headers in `cdhit/` name them. The algorithm, the data structures and the
output format are theirs. Questions about CD-HIT itself go to liwz@sdsc.edu or
https://github.com/weizhongli/cdhit, not here.

## Citation

Cite the original papers, not this port. These are the three the C++ prints in its `-h` text
(`cd_hit_ref1`, `cd_hit_ref2`, `cd_hit_ref3` in `src/cdhit_utility.rs`), and they are
reproduced there verbatim:

1. Weizhong Li & Adam Godzik. CD-HIT: a fast program for clustering and comparing large sets
   of protein or nucleotide sequences. Bioinformatics (2006) 22:1658-1659.
2. Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. CD-HIT: accelerated for
   clustering the next generation sequencing data. Bioinformatics (2012) 28:3150-3152.
3. Beifang Niu, Limin Fu, Shulei Sun & Weizhong Li. Artificial and natural duplicates in
   pyrosequencing reads of metagenomic data. BMC Bioinformatics (2010) 11:187.

(1) is the original method; (2) is the 4.x version this is translated from; (3) is the
duplicate detection in pyrosequencing reads that `cd-hit-454` implements. `cd-hit`,
`cd-hit-2d`, `cd-hit-est` and `cd-hit-est-2d` print (1) and (2); `cd-hit-454` prints all
three; `cd-hit-div` prints none.
