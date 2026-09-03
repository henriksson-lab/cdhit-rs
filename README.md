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

Peak RSS is slightly lower than the C++ too (89 MB vs 93 MB on the `-n 10` case), and no
longer grows with `-T`.

Since those numbers were taken, three further changes landed: a software-pipelined prefetch
of the word-table lookup, and an AVX-512 counting kernel for the nucleotide diagonal test
which takes the worst case -- `cd-hit-est -c 0.90`, where the word filter degenerates and
almost all time goes into that one function -- from parity to **1.83x**. Both are exact;
see `NAMING.md`.

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
