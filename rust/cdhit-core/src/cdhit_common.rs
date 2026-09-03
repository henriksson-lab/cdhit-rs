// =============================================================================
// CD-HI/CD-HIT
//
// Cluster Database at High Identity Threshold
//
// Rust translation of `cdhit-common.h` + `cdhit-common.c++` from CD-HIT 4.8.1.
//
// Original program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//
// modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
// =============================================================================
//
// Layout mirrors the C++: everything declared in `cdhit-common.h` comes first, then the
// definitions from `cdhit-common.c++` in their original order. See `../../NAMING.md` for
// the name mapping and the list of deliberate structural changes.

// This module is a line-by-line translation. Several lints fire on constructs that are
// faithful transcriptions of the C++ and are deliberately left alone:
//
// - `mut_range_bound`: `for j in from..imax_diag { ... from += 1 }` mirrors
//   `for (j=from; j<imax_diag; j++)`, where `from` is a separate variable and the bound is
//   fixed. Same for `for ks in i..m { i = ks + 1; }` in `do_clustering`.
// - `int_plus_one` / `manual_clamp` / `needless_range_loop` / `collapsible_if`: the C++
//   spelling is kept so the two sources can be read side by side.
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::mut_range_bound)]
#![allow(clippy::int_plus_one)]
#![allow(clippy::manual_clamp)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::if_same_then_else)]
#![allow(clippy::needless_bool)]
#![allow(clippy::large_enum_variant)]
#![allow(clippy::four_forward_slashes)]

use std::io::Write;


// ===========================================================================
// cdhit-common.h
// ===========================================================================

pub const CDHIT_VERSION: &str = "4.8.1";

/// Upstream compiles this in as a hard cap on sequence length. Here it survives only as the
/// scaling factor `ScoreMatrix` applies to every substitution score and gap penalty, which
/// is what makes the small `extra_score` diagonal bonus in `local_band_align` a tie-break
/// rather than a real contribution. Buffer sizes are derived from the input instead.
pub const MAX_SEQ: i32 = 655360;

pub const MAX_AA: usize = 23;
pub const MAX_NA: usize = 6;
/// Protein alphabet size. The EST front-ends override this to 4 with `#undef`/`#define`
/// before calling `init_naa`; here the value is simply passed to [`Naa::new`].
pub const MAX_UAA: i32 = 21;
pub const MAX_DES: usize = 300000;
pub const MAX_LINE_SIZE: usize = 300000;
pub const MAX_SEG: usize = 50;
pub const MAX_BIN_SWAP: f64 = 2e9;
pub const MAX_TABLE_SIZE: usize = 50000000;
pub const CLOCK_TICKS: i64 = 100;

pub const FAILED_FUNC: i32 = 1;
pub const OK_FUNC: i32 = 0;

pub const IS_REP: i16 = 1;
pub const IS_REDUNDANT: i16 = 2;
pub const IS_PROCESSED: i16 = 16;
pub const IS_MINUS_STRAND: i16 = 32;

/// `#define MAX_TABLE_SEQ 4000000`
pub const MAX_TABLE_SEQ: usize = 4000000;

/// `enum { DP_BACK_NONE=0, DP_BACK_LEFT_TOP=1, DP_BACK_LEFT=2, DP_BACK_TOP=3 }`
pub const DP_BACK_NONE: i32 = 0;
pub const DP_BACK_LEFT_TOP: i32 = 1;
pub const DP_BACK_LEFT: i32 = 2;
pub const DP_BACK_TOP: i32 = 3;

/// `typedef UINT4 INTs` under `#define LONG_SEQ`: word counts and AAP offsets are 32-bit.
pub type INTs = u32;

/// `struct IndexCount`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IndexCount {
    pub index: i32,
    pub count: i32,
}

impl IndexCount {
    pub fn new(index: i32, count: i32) -> Self {
        IndexCount { index, count }
    }
}

/// `class ScoreMatrix`
#[derive(Clone)]
pub struct ScoreMatrix {
    pub matrix: [[i32; MAX_AA]; MAX_AA],
    pub gap: i32,
    pub ext_gap: i32,
}

/// The file-scope `aa2idx` / `na2idx` arrays, which `setaa_to_na` and `-mask` mutate in
/// place. Held by value so a caller can run two clusterings with different alphabets.
#[derive(Clone, Copy)]
pub struct Alphabet {
    /// `aa2idx`: `letter - b'A'` -> residue index.
    pub aa2idx: [i32; 26],
    /// `na2idx`: `letter - b'A'` -> nucleotide index.
    pub na2idx: [i32; 26],
}

/// The `NAA1`..`NAA12` / `NAAN_array` globals that `InitNAA` fills in.
#[derive(Clone, Copy)]
pub struct Naa {
    /// `NAAN_array[k]` is the alphabet size raised to the `k`th power.
    pub naan_array: [i32; 13],
}

/// Capacity of a word's smallest block. With its one-slot header that makes the block 64
/// bytes, one cache line, so a word with up to seven entries costs a single miss.
const WORD_BLOCK0: u32 = 7;

/// Sentinel for "this word has no entries".
const NO_BLOCK: u32 = u32::MAX;

/// Capacity of a block of size class `k`.
#[inline]
fn word_block_cap(k: u32) -> u32 {
    WORD_BLOCK0 << k
}

/// Size class of a block with capacity `cap`.
#[inline]
fn word_block_class(cap: u32) -> u32 {
    (cap / WORD_BLOCK0).trailing_zeros()
}

/// `class WordTable`
pub struct WordTable {
    /// `indexCounts`: for each word code, the arena offset of its block of entries
    /// (`NO_BLOCK` when the word has none).
    ///
    /// Upstream gives every word its own `NVector<IndexCount>` -- a 16-byte header plus a
    /// separately `realloc`ed block. At protein's 4,084,101 words that is 65 MB of headers,
    /// all of it hit at random by `count_words`, the hottest loop in the program; a Rust
    /// `Vec` header would make it 98 MB. Keeping one 4-byte arena offset per word shrinks
    /// that random-access array 4x and removes the per-word allocation entirely, while a
    /// row's entries stay contiguous -- which matters just as much, because the words that
    /// dominate the loop are the ones with long entry lists.
    row_head: Vec<u32>,
    /// The block arena. A block at offset `o` is `[{ index: cap, count: len }, entries..]`,
    /// so a row's length is on the same cache line as its first entries.
    arena: Vec<IndexCount>,
    /// Blocks abandoned when a row outgrew them, by size class. A row that grows into class
    /// `k` takes one of these before extending the arena, so the doubling costs no memory
    /// beyond the largest block each row actually needs.
    free: Vec<Vec<u32>>,
    /// `sequences`: indices into `SequenceDB::sequences` (the C++ stores raw pointers).
    pub sequences: Vec<usize>,
    /// `NAA`: word length.
    pub naa: i32,
    /// `NAAN`: number of rows in the table.
    pub naan: i32,
    /// `is_aa`: non-zero for protein.
    pub is_aa: i8,
    pub size: usize,
    pub frag_count: i32,
}

/// `struct Options`
#[derive(Clone)]
pub struct Options {
    pub naa: i32,
    pub naan: i32,
    pub naa_top_limit: i32,

    /// `-M`, in bytes.
    pub max_memory: usize,
    /// `-l`
    pub min_length: i32,
    /// `-g`: false = the first match, true = the best match.
    pub cluster_best: bool,
    /// `-G`
    pub global_identity: bool,
    /// `-B`
    pub store_disk: bool,
    /// `-b`
    pub band_width: i32,
    /// `-c`
    pub cluster_thd: f64,
    /// `-D`
    pub distance_thd: f64,
    /// `-s`
    pub diff_cutoff: f64,
    /// `-s2`
    pub diff_cutoff2: f64,
    /// `-S`
    pub diff_cutoff_aa: i32,
    /// `-S2`
    pub diff_cutoff_aa2: i32,
    /// `-t`
    pub tolerance: i32,
    /// `-aL`
    pub long_coverage: f64,
    /// `-AL`
    pub long_control: i32,
    /// `-aS`
    pub short_coverage: f64,
    /// `-AS`
    pub short_control: i32,
    /// `-A`
    pub min_control: i32,
    /// `-uL`
    pub long_unmatch_per: f64,
    /// `-uS`
    pub short_unmatch_per: f64,
    /// `-U`
    pub unmatch_len: i32,
    /// `-D` in 454 mode.
    pub max_indel: i32,
    /// `-p`
    pub print: i32,
    /// `-d`
    pub des_len: i32,
    /// `-F`
    pub frag_size: i32,
    /// `-r`
    pub option_r: i32,
    /// `-T`
    pub threads: i32,
    /// `-P`
    pub pe_mode: i32,
    /// `-cx`
    pub trim_len: i32,
    /// `-cy`
    pub trim_len_r2: i32,
    /// `-ap`
    pub align_pos: i32,

    pub max_entries: usize,
    pub max_sequences: usize,
    pub mem_limit: usize,

    pub has_2d: bool,
    pub is_est: bool,
    pub is_454: bool,
    pub use_identity: bool,
    pub use_distance: bool,
    /// `-bak`
    pub backup_file: bool,

    pub input: String,
    pub input_pe: String,
    pub input2: String,
    pub input2_pe: String,
    pub output: String,
    pub output_pe: String,

    /// `-sc`
    pub sort_output: i32,
    /// `-sf`
    pub sort_outputf: i32,

    /// `-tmp`; the C++ keeps this in the file-scope `temp_dir`.
    pub temp_dir: String,
}

/// `struct Sequence`
pub struct Sequence {
    /// The residue-encoded sequence. Empty once `clear`ed or swapped out to disk.
    pub data: Vec<u8>,
    /// `size`: length of the sequence.
    pub size: i32,
    /// `size_R2`: for back-to-back merged paired-end reads, `size == size_R1 + size_R2`.
    pub size_r2: i32,

    /// `swap`: when set, the sequence lives in this temporary file rather than in `data`.
    pub swap: Option<usize>,
    /// Stream offset of the sequence within the swap file.
    pub offset: u64,

    /// Stream offset of the description string in the source database.
    pub des_begin: u64,
    pub des_begin2: u64,
    /// Total record length in the source database.
    pub tot_length: i32,
    pub tot_length2: i32,

    /// `identifier`: the full `>...` header line, leading `>` included.
    pub identifier: Vec<u8>,

    /// Index of the sequence in the original database.
    pub index: i32,
    pub state: i16,
    pub cluster_id: i32,
    pub identity: f32,
    pub distance: f32,
    pub coverage: [i32; 4],
}

/// `struct WorkingParam`
#[derive(Clone, Default)]
pub struct WorkingParam {
    pub aa1_cutoff: f64,
    /// `aas_cutoff` (a.k.a. aa2).
    pub aas_cutoff: f64,
    pub aan_cutoff: f64,
    pub len_upper_bound: i32,
    pub len_lower_bound: i32,

    pub len_eff: i32,
    pub aln_cover_flag: i32,
    /// `min_aln_lenS`
    pub min_aln_len_s: i32,
    /// `min_aln_lenL`
    pub min_aln_len_l: i32,
    pub required_aa1: i32,
    /// `required_aas` (a.k.a. aa2).
    pub required_aas: i32,
    pub required_aan: i32,
}

/// `lookCounts` plus its `size` field. `WordTable::count_words` resets the logical length to
/// zero without freeing, and relies on being able to write a terminating zero-count entry
/// one past the end, so the buffer keeps its own length separate from the `Vec`'s.
pub struct LookCounts {
    pub items: Vec<IndexCount>,
    pub size: usize,
}

/// `struct WorkingBuffer`
pub struct WorkingBuffer {
    pub taap: Vec<i32>,
    pub word_encodes: Vec<i32>,
    pub word_encodes_backup: Vec<i32>,
    pub word_encodes_no: Vec<INTs>,
    pub aap_list: Vec<INTs>,
    pub aap_begin: Vec<INTs>,
    /// `lookCounts`
    pub look_counts: LookCounts,
    /// `indexMapping`
    pub index_mapping: Vec<u32>,
    /// `score_mat` and `back_mat`, flattened to one allocation each with a uniform row
    /// stride. Upstream keeps a `Vector<NVector<...>>` and grows rows individually; a flat
    /// layout is equivalent because every cell `local_band_align` reads was written earlier
    /// in the same call (the band borders are initialised before the sweep, and each row
    /// only reads cells of itself and the row above that the sweep has already filled), so
    /// no state carries between calls. It also removes a pointer chase and a bounds check
    /// per cell access, and the per-row allocations.
    pub score_mat: Vec<i64>,
    pub back_mat: Vec<i32>,
    /// Row stride of `score_mat` / `back_mat`: the widest band seen so far.
    pub dp_stride: usize,
    /// `diag_score` and `diag_score2` packed into one word per diagonal: the low 32 bits
    /// are `diag_score[d]`, the high 32 are `diag_score2[d]`.
    ///
    /// Upstream keeps two `int` arrays, but every access touches the same index in both.
    /// Packing them halves the cache misses in `diag_test_*` -- where nearly all the time
    /// goes at low identity thresholds -- and turns the innermost update into a single
    /// 64-bit add instead of two read-modify-writes. Neither half can overflow into the
    /// other: `diag_score` counts word hits on one diagonal (at most `len2`), and
    /// `diag_score2` weights each by at most 4.
    pub diag_score: Vec<u64>,
    pub aan_list_comp: Vec<i32>,
    pub seqi_comp: Vec<u8>,
    pub total_bytes: i32,
}

/// `class SequenceDB`
#[derive(Default)]
pub struct SequenceDB {
    pub naan: i32,
    pub sequences: Vec<Sequence>,
    pub rep_seqs: Vec<i32>,

    pub total_letter: i64,
    pub total_desc: i64,
    pub max_len: usize,
    pub min_len: usize,
    pub len_n50: usize,
}

// ===========================================================================
// cdhit-common.c++
// ===========================================================================

/// `WITH_OPENMP`, the suffix stamped into the banner when OpenMP is compiled in.
pub fn with_openmp() -> &'static str {
    // The Rust port has no OpenMP build variant; multi-threading is unconditional.
    " (+OpenMP)"
}

/// `__DATE__`, in C's `"Mmm dd yyyy"` form with the day-of-month space-padded.
pub fn build_date() -> &'static str {
    env!("CDHIT_BUILD_DATE")
}

/// `__DATE__` with the day-padding space rewritten to `'0'`, which `SetOptions` does to
/// its local copy before printing the banner.
pub fn build_date_banner() -> String {
    let mut date: Vec<u8> = build_date().as_bytes().to_vec();
    for i in 1..date.len() {
        if date[i - 1] == b' ' && date[i] == b' ' {
            date[i] = b'0';
        }
    }
    String::from_utf8(date).unwrap()
}

/// `__TIME__`, in C's `"hh:mm:ss"` form.
pub fn build_time() -> &'static str {
    env!("CDHIT_BUILD_TIME")
}

/// `ctime(&tm)` for the current local time: `"Www Mmm dd hh:mm:ss yyyy\n"`, trailing
/// newline included.
pub fn ctime_now() -> String {
    chrono::Local::now()
        .format("%a %b %e %H:%M:%S %Y\n")
        .to_string()
}

/// `const char aa[]` -- the residue order the substitution matrices are laid out in.
pub const AA: &[u8; 23] = b"ARNDCQEGHILKMFPSTWYVBZX";

/// `aa2idx`, the default protein alphabet. `aa2idx[c - b'A']` is the index of residue `c`,
/// e.g. `aa2idx[b'A' - b'A'] == 0` and `aa2idx[b'M' - b'A'] == 12`.
///
/// ```text
/// idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
///          Q  R  S  T  U  V  W  X  Y  Z
/// ```
pub static AA2IDX: [i32; 26] = [
       0,    2,    4,    3,    6,   13,    7,    8,    9,   20,   11,   10,   12,
       2,   20,   14,    5,    1,   15,   16,   20,   19,   17,   20,   18,    6,
];

/// `BLOSUM62`, the lower triangle in row-major order as `ScoreMatrix::set_matrix` reads it.
///
/// ```text
/// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
/// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  2  6 20
/// ```
pub static BLOSUM62: [i32; 276] = [
       4,   -1,    5,   -2,    0,    6,   -2,   -2,    1,    6,    0,   -3,
      -3,   -3,    9,   -1,    1,    0,    0,   -3,    5,   -1,    0,    0,
       2,   -4,    2,    5,    0,   -2,    0,   -1,   -3,   -2,   -2,    6,
      -2,    0,    1,   -1,   -3,    0,    0,   -2,    8,   -1,   -3,   -3,
      -3,   -1,   -3,   -3,   -4,   -3,    4,   -1,   -2,   -3,   -4,   -1,
      -2,   -3,   -4,   -3,    2,    4,   -1,    2,    0,   -1,   -3,    1,
       1,   -2,   -1,   -3,   -2,    5,   -1,   -1,   -2,   -3,   -1,    0,
      -2,   -3,   -2,    1,    2,   -1,    5,   -2,   -3,   -3,   -3,   -2,
      -3,   -3,   -3,   -1,    0,    0,   -3,    0,    6,   -1,   -2,   -2,
      -1,   -3,   -1,   -1,   -2,   -2,   -3,   -3,   -1,   -2,   -4,    7,
       1,   -1,    1,    0,   -1,    0,    0,    0,   -1,   -2,   -2,    0,
      -1,   -2,   -1,    4,    0,   -1,    0,   -1,   -1,   -1,   -1,   -2,
      -2,   -1,   -1,   -1,   -1,   -2,   -1,    1,    5,   -3,   -3,   -4,
      -4,   -2,   -2,   -3,   -2,   -2,   -3,   -2,   -3,   -1,    1,   -4,
      -3,   -2,   11,   -2,   -2,   -2,   -3,   -2,   -1,   -2,   -3,    2,
      -1,   -1,   -2,   -1,    3,   -3,   -2,   -2,    2,    7,    0,   -3,
      -3,   -3,   -1,   -2,   -2,   -3,   -3,    3,    1,   -2,    1,   -1,
      -2,   -2,    0,   -3,   -1,    4,   -2,   -1,    3,    4,   -3,    0,
       1,   -1,    0,   -3,   -4,    0,   -3,   -3,   -2,    0,   -1,   -4,
      -3,   -3,    4,   -1,    0,    0,    1,   -3,    3,    4,   -2,    0,
      -3,   -3,    1,   -1,   -3,   -1,    0,   -1,   -3,   -2,   -2,    1,
       4,    0,   -1,   -1,   -1,   -2,   -1,   -1,   -1,   -1,   -1,   -1,
      -1,   -1,   -1,   -2,    0,    0,   -2,   -1,   -1,   -1,   -1,   -1,
];

/// `na2idx`, the default nucleotide alphabet: A,C,G,T/U -> 0..=3, everything else -> 4.
/// `-mask` rewrites the listed letters to 5.
pub static NA2IDX: [i32; 26] = [
       0,    4,    1,    4,    4,    4,    2,    4,    4,    4,    4,    4,    4,
       4,    4,    4,    4,    4,    4,    3,    3,    4,    4,    4,    4,    4,
];

/// `BLOSUM62_na`, the lower triangle of the nucleotide substitution matrix.
///
/// ```text
/// A  C  G  T  U  N  X
/// 0  1  2  3  3  4  5
/// ```
pub static BLOSUM62_NA: [i32; 28] = [
       2,   -2,    2,   -2,   -2,    2,   -2,
      -2,   -2,    2,   -2,   -2,   -2,    1,
       2,   -2,   -2,   -2,   -2,   -2,    1,
       0,    0,    0,    0,    0,    0,    1,
];

// ---------------------------------------------------------------------------
// TempFile / TempFiles / OpenTempFile / CleanUpTempFiles
//
// COMPAT: dead code upstream. `OpenTempFile` has no callers in 4.8.1 and nothing ever
// assigns `Sequence::swap`, so the registry is always empty: `bomb_error`'s
// `temp_files.Clear()` and the `atexit(CleanUpTempFiles)` handler are both no-ops, and
// `-B` (`store_disk`) does nothing. `Sequence::swap_in` / `swap_out` are translated below
// and are equally inert. Not reproduced here.
// ---------------------------------------------------------------------------

impl Naa {
    /// `InitNAA`: fill `NAAN_array[k] = max^k`.
    ///
    /// COMPAT: upstream computes these as `int`, so for `max = 21` everything from
    /// `NAA8` up overflows. That is harmless only because the protein front-ends cap the
    /// word length at 5 and the EST front-ends `#undef MAX_UAA` to 4 before calling this.
    /// The wrapping arithmetic is preserved so any table size matches upstream exactly.
    pub fn new(max: i32) -> Naa {
        let mut a = [0i32; 13];
        a[0] = 1;
        a[1] = max;
        a[2] = a[1].wrapping_mul(a[1]);
        a[3] = a[1].wrapping_mul(a[2]);
        a[4] = a[2].wrapping_mul(a[2]);
        a[5] = a[2].wrapping_mul(a[3]);
        a[6] = a[3].wrapping_mul(a[3]);
        a[7] = a[3].wrapping_mul(a[4]);
        a[8] = a[4].wrapping_mul(a[4]);
        a[9] = a[4].wrapping_mul(a[5]);
        a[10] = a[5].wrapping_mul(a[5]);
        a[11] = a[5].wrapping_mul(a[6]);
        a[12] = a[6].wrapping_mul(a[6]);
        Naa { naan_array: a }
    }

    /// `NAA1`, the alphabet size.
    pub fn naa1(&self) -> i32 {
        self.naan_array[1]
    }
    /// `NAA2`
    pub fn naa2(&self) -> i32 {
        self.naan_array[2]
    }
    /// `NAA3`
    pub fn naa3(&self) -> i32 {
        self.naan_array[3]
    }
    /// `NAA4`
    pub fn naa4(&self) -> i32 {
        self.naan_array[4]
    }
}

/// `make_comp_iseq`: write the reverse complement of `iseq` into `iseq_comp`.
pub fn make_comp_iseq(len: usize, iseq_comp: &mut [u8], iseq: &[u8]) {
    let c: [u8; 6] = [3, 2, 1, 0, 4, 5];
    for i in 0..len {
        iseq_comp[i] = c[iseq[len - i - 1] as usize];
    }
}

impl Default for Options {
    /// `Options::Options()`
    fn default() -> Options {
        Options {
            naa: 5,
            naan: 0,
            naa_top_limit: 5,
            max_memory: 800000000,
            min_length: 10,
            cluster_best: false,
            global_identity: true,
            store_disk: false,
            band_width: 20,
            cluster_thd: 0.9,
            distance_thd: 0.0,
            diff_cutoff: 0.0,
            diff_cutoff2: 1.0,
            diff_cutoff_aa: 99999999,
            diff_cutoff_aa2: 0,
            tolerance: 2,
            long_coverage: 0.0,
            long_control: 99999999,
            short_coverage: 0.0,
            short_control: 99999999,
            min_control: 0,
            long_unmatch_per: 1.0,
            short_unmatch_per: 1.0,
            unmatch_len: 99999999,
            max_indel: 1,
            print: 0,
            des_len: 20,
            frag_size: 0,
            option_r: 1,
            threads: 1,
            pe_mode: 0,
            trim_len: 0,
            trim_len_r2: 0,
            align_pos: 0,
            max_entries: 0,
            max_sequences: 1 << 20,
            mem_limit: 100000000,
            has_2d: false,
            is_est: false,
            is_454: false,
            use_identity: false,
            use_distance: false,
            backup_file: false,
            input: String::new(),
            input_pe: String::new(),
            input2: String::new(),
            input2_pe: String::new(),
            output: String::new(),
            output_pe: String::new(),
            sort_output: 0,
            sort_outputf: 0,
            temp_dir: String::new(),
        }
    }
}

/// `atoi`: C semantics -- leading whitespace, optional sign, digits, and 0 on anything else.
fn atoi(value: &str) -> i32 {
    let s = value.trim_start();
    let (sign, rest) = match s.strip_prefix('-') {
        Some(r) => (-1i64, r),
        None => (1i64, s.strip_prefix('+').unwrap_or(s)),
    };
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
    let n: i64 = digits.parse().unwrap_or(0);
    (sign * n).clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

/// `atoll`
fn atoll(value: &str) -> i64 {
    let s = value.trim_start();
    let (sign, rest) = match s.strip_prefix('-') {
        Some(r) => (-1i128, r),
        None => (1i128, s.strip_prefix('+').unwrap_or(s)),
    };
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
    let n: i128 = digits.parse().unwrap_or(0);
    (sign * n).clamp(i64::MIN as i128, i64::MAX as i128) as i64
}

/// `atof`
fn atof(value: &str) -> f64 {
    let s = value.trim_start();
    let mut end = 0;
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let mut seen_digit = false;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        seen_digit = true;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            seen_digit = true;
        }
    }
    if seen_digit {
        end = i;
        if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
            let mut j = i + 1;
            if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
                j += 1;
            }
            let mut exp_digit = false;
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
                exp_digit = true;
            }
            if exp_digit {
                end = j;
            }
        }
    }
    s[..end].parse().unwrap_or(0.0)
}

impl Options {
    /// `Options::SetOptionCommon`
    pub fn set_option_common(&mut self, flag: &str, value: &str, num_procs: i32) -> bool {
        let intval = atoi(value);
        match flag {
            "-i" => self.input = value.to_string(),
            "-j" => self.input_pe = value.to_string(),
            "-o" => self.output = value.to_string(),
            "-op" => self.output_pe = value.to_string(),
            "-M" => self.max_memory = (atoll(value) * 1000000) as usize,
            "-l" => self.min_length = intval,
            "-c" => {
                self.cluster_thd = atof(value);
                self.use_identity = true;
            }
            "-D" => {
                self.distance_thd = atof(value);
                self.use_distance = true;
            }
            "-b" => self.band_width = intval,
            "-n" => self.naa = intval,
            "-d" => self.des_len = intval,
            "-s" => self.diff_cutoff = atof(value),
            "-S" => self.diff_cutoff_aa = intval,
            "-B" => self.store_disk = intval != 0,
            "-P" => self.pe_mode = intval,
            "-cx" => self.trim_len = intval,
            "-cy" => self.trim_len_r2 = intval,
            "-ap" => self.align_pos = intval,
            "-sc" => self.sort_output = intval,
            "-sf" => self.sort_outputf = intval,
            "-p" => self.print = intval,
            "-g" => self.cluster_best = intval != 0,
            "-G" => self.global_identity = intval != 0,
            "-aL" => self.long_coverage = atof(value),
            "-AL" => self.long_control = intval,
            "-aS" => self.short_coverage = atof(value),
            "-AS" => self.short_control = intval,
            "-A" => self.min_control = intval,
            "-uL" => self.long_unmatch_per = atof(value),
            "-uS" => self.short_unmatch_per = atof(value),
            "-U" => self.unmatch_len = intval,
            "-tmp" => self.temp_dir = value.to_string(),
            "-bak" => self.backup_file = intval != 0,
            "-T" => {
                let cpu = num_procs;
                self.threads = intval;
                if self.threads > cpu {
                    self.threads = cpu;
                    println!("Warning: total number of CPUs in the system is {}", cpu);
                } else if self.threads < 0 {
                    self.threads += cpu;
                    if self.threads < 0 {
                        self.threads = 0;
                    }
                }
                if self.threads == 0 {
                    self.threads = cpu;
                    println!("total number of CPUs in the system is {}", cpu);
                }
                if self.threads != intval {
                    println!("Actual number of CPUs to be used: {}\n", self.threads);
                }
            }
            _ => return false,
        }
        true
    }

    /// `Options::SetOption`
    pub fn set_option(
        &mut self,
        flag: &str,
        value: &str,
        num_procs: i32,
        mat: &mut ScoreMatrix,
        alphabet: &mut Alphabet,
    ) -> bool {
        if self.is_454 {
            match flag {
                "-s" | "-S" | "-G" | "-A" | "-r" => return false,
                "-D" => {
                    self.max_indel = atoi(value);
                    return true;
                }
                _ => {}
            }
        }
        if self.set_option_common(flag, value, num_procs) {
            return true;
        }
        if flag == "-t" {
            self.tolerance = atoi(value);
        } else if flag == "-F" {
            self.frag_size = atoi(value);
        } else if self.has_2d && self.set_option_2d(flag, value, num_procs) {
            return true;
        } else if self.is_est && self.set_option_est(flag, value, num_procs, mat, alphabet) {
            return true;
        } else {
            return false;
        }
        true
    }

    /// `Options::SetOption2D`
    pub fn set_option_2d(&mut self, flag: &str, value: &str, num_procs: i32) -> bool {
        if self.set_option_common(flag, value, num_procs) {
            return true;
        }
        match flag {
            "-i2" => self.input2 = value.to_string(),
            "-j2" => self.input2_pe = value.to_string(),
            "-s2" => self.diff_cutoff2 = atof(value),
            "-S2" => self.diff_cutoff_aa2 = atoi(value),
            _ => return false,
        }
        true
    }

    /// `Options::SetOptionEST`
    pub fn set_option_est(
        &mut self,
        flag: &str,
        value: &str,
        num_procs: i32,
        mat: &mut ScoreMatrix,
        alphabet: &mut Alphabet,
    ) -> bool {
        self.naa_top_limit = 12;
        if self.set_option_common(flag, value, num_procs) {
            return true;
        }
        match flag {
            "-r" => self.option_r = atoi(value),
            "-gap" => mat.gap = MAX_SEQ.wrapping_mul(atoi(value)),
            "-gap-ext" => mat.ext_gap = MAX_SEQ.wrapping_mul(atoi(value)),
            "-match" => mat.set_match(atoi(value)),
            "-mismatch" => mat.set_mismatch(atoi(value)),
            "-mask" => {
                for ch in value.bytes() {
                    let ch = ch.to_ascii_uppercase();
                    if !ch.is_ascii_uppercase() {
                        continue;
                    }
                    alphabet.na2idx[(ch - b'A') as usize] = 5;
                }
                alphabet.setaa_to_na();
            }
            _ => return false,
        }
        true
    }

    /// `Options::SetOptions`
    ///
    /// `argv` is the full argument vector, program name included, as in C.
    pub fn set_options(
        &mut self,
        argv: &[String],
        twod: bool,
        est: bool,
        num_procs: i32,
        mat: &mut ScoreMatrix,
        alphabet: &mut Alphabet,
    ) -> bool {
        let argc = argv.len();
        println!("================================================================");
        println!(
            "Program: CD-HIT, V{}{}, {}, {}",
            CDHIT_VERSION,
            with_openmp(),
            build_date_banner(),
            build_time()
        );
        print!("Command:");
        let mut n = 9;
        for a in argv.iter() {
            n += a.len() + 1;
            if n >= 64 {
                print!("\n         {}", a);
                n = a.len() + 9;
            } else {
                print!(" {}", a);
            }
        }
        println!("\n");
        print!("Started: {}", ctime_now());
        println!("================================================================");
        println!("                            Output                              ");
        println!("----------------------------------------------------------------");
        self.has_2d = twod;
        self.is_est = est;
        let mut i = 1;
        while i + 1 < argc {
            if !self.set_option(&argv[i], &argv[i + 1], num_procs, mat, alphabet) {
                return false;
            }
            i += 2;
        }
        if i < argc {
            return false;
        }
        // COMPAT: `atexit(CleanUpTempFiles)` -- inert, see the TempFile note above.
        true
    }
}

impl Options {
    /// `Options::Validate`
    pub fn validate(&mut self) {
        if self.use_identity && self.use_distance {
            bomb_error("can not use both identity cutoff and distance cutoff");
        }
        if self.use_distance {
            if (self.distance_thd > 1.0) || (self.distance_thd < 0.0) {
                bomb_error("invalid distance threshold");
            }
        } else if self.is_est {
            if (self.cluster_thd > 1.0) || (self.cluster_thd < 0.8) {
                bomb_error("invalid clstr threshold, should >=0.8");
            }
        } else if (self.cluster_thd > 1.0) || (self.cluster_thd < 0.4) {
            bomb_error("invalid clstr");
        }

        if self.input.is_empty() {
            bomb_error("no input file");
        }
        if self.output.is_empty() {
            bomb_error("no output file");
        }
        if self.pe_mode != 0 {
            if self.input_pe.is_empty() {
                bomb_error("no input file for R2 sequences in PE mode");
            }
            if self.output_pe.is_empty() {
                bomb_error("no output file for R2 sequences in PE mode");
            }
        }
        if self.is_est && self.align_pos == 1 {
            self.option_r = 0;
        }

        if self.band_width < 1 {
            bomb_error("invalid band width");
        }
        if self.naa < 2 || self.naa > self.naa_top_limit {
            bomb_error("invalid word length");
        }
        if self.des_len < 0 {
            bomb_error("too short description, not enough to identify sequences");
        }
        if !self.is_est && (self.tolerance < 0 || self.tolerance > 5) {
            bomb_error("invalid tolerance");
        }
        if (self.diff_cutoff < 0.0) || (self.diff_cutoff > 1.0) {
            bomb_error("invalid value for -s");
        }
        if self.diff_cutoff_aa < 0 {
            bomb_error("invalid value for -S");
        }
        if self.has_2d {
            if (self.diff_cutoff2 < 0.0) || (self.diff_cutoff2 > 1.0) {
                bomb_error("invalid value for -s2");
            }
            if self.diff_cutoff_aa2 < 0 {
                bomb_error("invalid value for -S2");
            }
            if self.pe_mode != 0 && self.input2_pe.is_empty() {
                bomb_error("no input file for R2 sequences for 2nd db in PE mode");
            }
        }
        if !self.global_identity {
            self.print = 1;
        }
        if self.short_coverage < self.long_coverage {
            self.short_coverage = self.long_coverage;
        }
        if self.short_control > self.long_control {
            self.short_control = self.long_control;
        }
        if !self.global_identity && (self.short_coverage == 0.0) && (self.min_control == 0) {
            bomb_error("You are using local identity, but no -aS -aL -A option");
        }
        if self.frag_size < 0 {
            bomb_error("invalid fragment size");
        }

        let message = "Your word length is {}, using {} may be faster!";
        let _ = message;
        if !self.is_est && self.tolerance != 0 {
            let clstr_idx = (self.cluster_thd * 100.0) as i32 - NAA_STAT_START_PERCENT;
            let tcutoff =
                NAA_STAT[(self.tolerance - 1) as usize][clstr_idx as usize][(5 - self.naa) as usize];

            if tcutoff < 5 {
                bomb_error(
                    "Too low cluster threshold for the word length.\n\
                     Increase the threshold or the tolerance, or decrease the word length.",
                );
            }
            let mut i = 5;
            while i > self.naa {
                if NAA_STAT[(self.tolerance - 1) as usize][clstr_idx as usize][(5 - i) as usize] > 10
                {
                    println!("Your word length is {}, using {} may be faster!", self.naa, i);
                    break;
                }
                i -= 1;
            }
        } else if self.is_est {
            if self.cluster_thd > 0.9 && self.naa < 8 {
                println!("Your word length is {}, using {} may be faster!", self.naa, 8);
            } else if self.cluster_thd > 0.87 && self.naa < 5 {
                println!("Your word length is {}, using {} may be faster!", self.naa, 5);
            } else if self.cluster_thd > 0.80 && self.naa < 4 {
                println!("Your word length is {}, using {} may be faster!", self.naa, 4);
            } else if self.cluster_thd > 0.75 && self.naa < 3 {
                println!("Your word length is {}, using {} may be faster!", self.naa, 3);
            }
        } else if self.cluster_thd > 0.85 && self.naa < 5 {
            println!("Your word length is {}, using {} may be faster!", self.naa, 5);
        } else if self.cluster_thd > 0.80 && self.naa < 4 {
            println!("Your word length is {}, using {} may be faster!", self.naa, 4);
        } else if self.cluster_thd > 0.75 && self.naa < 3 {
            println!("Your word length is {}, using {} may be faster!", self.naa, 3);
        }

        if (self.min_length + 1) < self.naa {
            bomb_error("Too short -l, redefine it");
        }
    }

    /// `Options::Print`
    pub fn print(&self) {
        println!("isEST = {}", self.is_est as i32);
        println!("has2D = {}", self.has_2d as i32);
        println!("NAA = {}", self.naa);
        println!("NAA_top_limit = {}", self.naa_top_limit);
        println!("min_length = {}", self.min_length);
        println!("cluster_best = {}", self.cluster_best as i32);
        println!("global_identity = {}", self.global_identity as i32);
        println!("cluster_thd = {}", self.cluster_thd);
        println!("diff_cutoff = {}", self.diff_cutoff);
        println!("diff_cutoff_aa = {}", self.diff_cutoff_aa);
        println!("tolerance = {}", self.tolerance);
        println!("long_coverage = {}", self.long_coverage);
        println!("long_control = {}", self.long_control);
        println!("short_coverage = {}", self.short_coverage);
        println!("short_control = {}", self.short_control);
        println!("frag_size = {}", self.frag_size);
        println!("option_r = {}", self.option_r);
        println!("print = {}", self.print);
    }
}

/// `bomb_error`
pub fn bomb_error(message: &str) -> ! {
    eprintln!("\nFatal Error:\n{}\nProgram halted !!\n", message);
    std::process::exit(1);
}

/// `bomb_error( const char*, const char* )`
pub fn bomb_error2(message: &str, message2: &str) -> ! {
    eprintln!(
        "\nFatal Error:\n{} {}\nProgram halted !!\n",
        message, message2
    );
    std::process::exit(1);
}

/// `bomb_warning`
pub fn bomb_warning(message: &str) {
    eprintln!("\nWarning:\n{}\nNot fatal, but may affect results !!\n", message);
}

/// `bomb_warning( const char*, const char* )`
pub fn bomb_warning2(message: &str, message2: &str) {
    eprintln!(
        "\nWarning:\n{} {}\nNot fatal, but may affect results !!\n",
        message, message2
    );
}

/// C's `isspace` in the "C" locale: space, `\t`, `\n`, `\v`, `\f`, `\r`. Rust's
/// `u8::is_ascii_whitespace` omits `\v`, so it cannot be used here.
#[inline]
pub fn c_isspace(ch: u8) -> bool {
    matches!(ch, b' ' | 0x09..=0x0d)
}

/// `format_seq`: upper-case in place and drop everything that is not a letter.
pub fn format_seq(seq: &mut Vec<u8>) {
    let mut j = 0;
    for i in 0..seq.len() {
        let c1 = seq[i].to_ascii_uppercase();
        if c1.is_ascii_alphabetic() {
            seq[j] = c1;
            j += 1;
        }
    }
    seq.truncate(j);
}

/// `strrev`
pub fn strrev(p: &mut [u8]) {
    p.reverse();
}

/// Split a packed diagonal cell into `(diag_score, diag_score2)`.
#[inline]
fn unpack_diag(v: u64) -> (i32, i32) {
    (v as u32 as i32, (v >> 32) as i32)
}

////For smiple len1 <= len2, len2 is for existing representative
////walk along all diag path of two sequences,
////find the diags with most aap
////return top n diags
////added on 2006 11 13
////band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                            XXXXXXXXXXXXXXX                  seq1
////band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                             XXXXXXXXXXXXXXX                 seq1
////extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
////    band = len2-1                            XXXXXXXXXXXXXXX seq1
////band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                           XXXXXXXXXXXXXXX                   seq1
////extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
////              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
////index of diag_score = band+len1-1;
/// `diag_test_aapn`
///
/// `buffer.taap` / `aap_begin` / `aap_list` must already describe seq1's dipeptides, i.e.
/// [`WorkingBuffer::compute_aap`] must have been called for the query.
pub fn diag_test_aapn(
    naa1: i32,
    iseq2: &[u8],
    len1: i32,
    len2: i32,
    buffer: &mut WorkingBuffer,
    best_sum: &mut i32,
    band_width: i32,
    band_left: &mut i32,
    band_center: &mut i32,
    band_right: &mut i32,
    required_aa1: i32,
    options: &Options,
) -> i32 {
    let nall = len1 + len2 - 1;
    let taap = &buffer.taap;
    let aap_begin = &buffer.aap_begin;
    let aap_list = &buffer.aap_list;
    let diag_score = &mut buffer.diag_score;

    if nall as usize > diag_score.len() {
        bomb_error("in diag_test_aapn, MAX_DIAG reached");
    }

    //find the best band range
    //  int band_b = required_aa1;
    let band_b = if required_aa1 - 1 >= 0 { required_aa1 - 1 } else { 0 }; // on dec 21 2001
    let band_e = nall - band_b;

    // Upstream clears all `nall` diagonals before counting. Only `[band_b, band_e]` is ever
    // read back, though -- the scans below never look outside it -- so clearing just that
    // range is equivalent and, at low identity thresholds where `band_b` is large and every
    // representative becomes a candidate, is most of the work this function does. Counts
    // that land outside the range are written and never read; the next call clears whatever
    // part of the array it is going to read.
    if band_b <= band_e {
        diag_score[band_b as usize..=band_e as usize].fill(0);
    }

    let mut i1 = len1 - 1;
    for w in iseq2[..len2 as usize].windows(2) {
        let c22 = ((w[0] as i32) * naa1 + w[1] as i32) as usize;
        debug_assert!(c22 < taap.len() && c22 < aap_begin.len());
        // SAFETY: residues are < naa1, so `c22 < naa1^2`, and `WorkingBuffer::set` sizes
        // both arrays to at least `MAX_UAA^2` (441).
        let j = unsafe { *taap.get_unchecked(c22) };
        if j == 0 {
            i1 += 1;
            continue;
        }
        let cpx = 1 + (w[0] != w[1]) as i32;
        let step = 1u64 | ((cpx as u64) << 32);
        // SAFETY: as above.
        let m = unsafe { *aap_begin.get_unchecked(c22) } as usize;
        // SAFETY: `compute_aap`/`compute_aap2` fill `aap_list` so that
        // `aap_begin[c22] + taap[c22]` never exceeds the number of words in seq1, which is
        // at most `aap_list.len()`.
        debug_assert!(m + j as usize <= aap_list.len());
        for &b in unsafe { aap_list.get_unchecked(m..m + j as usize) } {
            let d = (i1 - b as i32) as usize;
            debug_assert!(d < nall as usize);
            // SAFETY: `aap_list` holds positions in seq1, so `0 <= b <= len1 - 2`, and
            // `len1 - 1 <= i1 <= len1 + len2 - 3`. Hence `1 <= d <= len1 + len2 - 3 < nall`,
            // and the array was just sized and zeroed to `nall`. This is the innermost loop
            // of the whole program at low identity thresholds.
            unsafe { *diag_score.get_unchecked_mut(d) += step };
        }
        i1 += 1;
    }

    let band_m = if band_b + band_width - 1 < band_e {
        band_b + band_width - 1
    } else {
        band_e
    };
    let mut best_score = 0;
    let mut best_score2 = 0;
    let mut max_diag2 = 0;
    let mut imax_diag = 0;
    // Every index below is at most `band_e`, and `band_e <= nall`, so all of them are in
    // the range that was just zeroed. (`band_m` equals `band_e` only when `band_b > 0`,
    // which puts `band_e` strictly below `nall`.) The sliding scan runs once per diagonal
    // for every candidate pair and is most of the program's time at low identity
    // thresholds, hence the unchecked access.
    debug_assert!((band_e as usize) < diag_score.len());
    let ds = &diag_score[..band_e as usize + 1];
    for i in band_b..=band_m {
        // SAFETY: `band_m <= band_e` and `ds` has `band_e + 1` elements.
        let e = unsafe { *ds.get_unchecked(i as usize) };
        let (e0, e1) = unpack_diag(e);
        best_score += e0;
        best_score2 += e1;
        if e1 > max_diag2 {
            max_diag2 = e1;
            imax_diag = i;
        }
    }
    let mut from = band_b;
    let mut end = band_m;
    let mut score = best_score;
    let mut score2 = best_score2;
    let mut k = from;
    let mut j = band_m + 1;
    while j < band_e {
        // SAFETY: `k < j < band_e` and `ds` has `band_e + 1` elements.
        let (k0, k1) = unpack_diag(unsafe { *ds.get_unchecked(k as usize) });
        let (j0, j1) = unpack_diag(unsafe { *ds.get_unchecked(j as usize) });
        score -= k0;
        score += j0;
        score2 -= k1;
        score2 += j1;
        if score2 > best_score2 {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if j1 > max_diag2 {
                max_diag2 = j1;
                imax_diag = j;
            }
        }
        j += 1;
        k += 1;
    }
    let mlen = if imax_diag > len1 { nall - imax_diag } else { imax_diag };
    let emax = ((1.0 - options.cluster_thd) * mlen as f64) as i32 + 1;
    for j in from..imax_diag {
        // if aap pairs fail to open gap
        if (imax_diag - j) > emax || unpack_diag(ds[j as usize]).0 < 1 {
            best_score -= unpack_diag(ds[j as usize]).0;
            from += 1;
        } else {
            break;
        }
    }
    let mut j = end;
    while j > imax_diag {
        // if aap pairs fail to open gap
        if (j - imax_diag) > emax || unpack_diag(ds[j as usize]).0 < 1 {
            best_score -= unpack_diag(ds[j as usize]).0;
            end -= 1;
        } else {
            break;
        }
        j -= 1;
    }

    *band_left = from - len1 + 1;
    *band_right = end - len1 + 1;
    *band_center = imax_diag - len1 + 1;
    *best_sum = best_score;
    OK_FUNC
}
// END diag_test_aapn

/// `diag_test_aapn_est`
///
/// `buffer.taap` / `aap_begin` / `aap_list` must already describe seq1's 4-mers, i.e.
/// [`WorkingBuffer::compute_aap2`] must have been called for the query.
pub fn diag_test_aapn_est(
    naa1: i32,
    iseq2: &[u8],
    len1: i32,
    len2: i32,
    buffer: &mut WorkingBuffer,
    best_sum: &mut i32,
    band_width: i32,
    band_left: &mut i32,
    band_center: &mut i32,
    band_right: &mut i32,
    required_aa1: i32,
    options: &Options,
) -> i32 {
    let nall = len1 + len2 - 1;
    let naa2 = naa1 * naa1;
    let naa3 = naa2 * naa1;
    let taap = &buffer.taap;
    let aap_begin = &buffer.aap_begin;
    let aap_list = &buffer.aap_list;
    let diag_score = &mut buffer.diag_score;

    if nall as usize > diag_score.len() {
        bomb_error("in diag_test_aapn_est, MAX_DIAG reached");
    }

    //find the best band range
    //  int band_b = required_aa1;
    let mut band_b = if required_aa1 - 1 >= 0 { required_aa1 - 1 } else { 0 }; // on dec 21 2001
    let mut band_e = nall - band_b;

    if options.is_454 {
        band_b = len1 - band_width;
        band_e = len1 + band_width;
        if band_b < 0 {
            band_b = 0;
        }
        if band_e > nall {
            band_e = nall;
        }
    }

    // Upstream clears all `nall` diagonals before counting. Only `[band_b, band_e]` is ever
    // read back, though -- the scans below never look outside it -- so clearing just that
    // range is equivalent and, at low identity thresholds where `band_b` is large and every
    // representative becomes a candidate, is most of the work this function does. Counts
    // that land outside the range are written and never read; the next call clears whatever
    // part of the array it is going to read.
    if band_b <= band_e {
        diag_score[band_b as usize..=band_e as usize].fill(0);
    }

    let mut i1 = len1 - 1;
    for w in iseq2[..len2 as usize].windows(4) {
        // Read the 4-mer as one 32-bit word. Bases are 0..=3, N is 4 and the mask code is
        // 5, so a byte to skip is exactly one with a bit above the low two set: one mask
        // and one test replaces four comparisons.
        let word = u32::from_le_bytes([w[0], w[1], w[2], w[3]]);
        if word & 0xfcfc_fcfc != 0 {
            i1 += 1;
            continue; //skip N
        }
        let c0 = (word & 0xff) as i32;
        let c1 = ((word >> 8) & 0xff) as i32;
        let c2 = ((word >> 16) & 0xff) as i32;
        let c3 = (word >> 24) as i32;

        let c22 = (c0 * naa3 + c1 * naa2 + c2 * naa1 + c3) as usize;
        debug_assert!(c22 < taap.len() && c22 < aap_begin.len());
        // SAFETY: each base is < naa1, so `c22 < naa1^4`, and `WorkingBuffer::set` sizes
        // both arrays to `MAX_UAA^4` (194481), comfortably above the nucleotide 4^4 = 256.
        let j = unsafe { *taap.get_unchecked(c22) };
        if j == 0 {
            i1 += 1;
            continue;
        }
        let cpx = 1 + (c0 != c1) as i32 + (c1 != c2) as i32 + (c2 != c3) as i32;
        let step = 1u64 | ((cpx as u64) << 32);
        // SAFETY: as above.
        let m = unsafe { *aap_begin.get_unchecked(c22) } as usize;
        // SAFETY: `compute_aap`/`compute_aap2` fill `aap_list` so that
        // `aap_begin[c22] + taap[c22]` never exceeds the number of words in seq1, which is
        // at most `aap_list.len()`.
        debug_assert!(m + j as usize <= aap_list.len());
        for &b in unsafe { aap_list.get_unchecked(m..m + j as usize) } {
            let d = (i1 - b as i32) as usize;
            debug_assert!(d < nall as usize);
            // SAFETY: `aap_list` holds positions in seq1, so `0 <= b <= len1 - 4`, and
            // `len1 - 1 <= i1 <= len1 + len2 - 5`. Hence `3 <= d <= len1 + len2 - 5 < nall`,
            // and the array was just sized and zeroed to `nall`. Nearly all of the run time
            // at low identity thresholds is in this loop.
            unsafe { *diag_score.get_unchecked_mut(d) += step };
        }
        i1 += 1;
    }

    let band_m = if band_b + band_width - 1 < band_e {
        band_b + band_width - 1
    } else {
        band_e
    };
    let mut best_score = 0;
    let mut best_score2 = 0;
    let mut max_diag2 = 0;
    let mut imax_diag = 0;
    // Every index below is at most `band_e`, and `band_e <= nall`, so all of them are in
    // the range that was just zeroed. (`band_m` equals `band_e` only when `band_b > 0`,
    // which puts `band_e` strictly below `nall`.) The sliding scan runs once per diagonal
    // for every candidate pair and is most of the program's time at low identity
    // thresholds, hence the unchecked access.
    debug_assert!((band_e as usize) < diag_score.len());
    let ds = &diag_score[..band_e as usize + 1];
    for i in band_b..=band_m {
        // SAFETY: `band_m <= band_e` and `ds` has `band_e + 1` elements.
        let e = unsafe { *ds.get_unchecked(i as usize) };
        let (e0, e1) = unpack_diag(e);
        best_score += e0;
        best_score2 += e1;
        if e1 > max_diag2 {
            max_diag2 = e1;
            imax_diag = i;
        }
    }
    let mut from = band_b;
    let mut end = band_m;
    let mut score = best_score;
    let mut score2 = best_score2;
    let mut k = from;
    let mut j = band_m + 1;
    while j < band_e {
        // SAFETY: `k < j < band_e` and `ds` has `band_e + 1` elements.
        let (k0, k1) = unpack_diag(unsafe { *ds.get_unchecked(k as usize) });
        let (j0, j1) = unpack_diag(unsafe { *ds.get_unchecked(j as usize) });
        score -= k0;
        score += j0;
        score2 -= k1;
        score2 += j1;
        if score2 > best_score2 {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if j1 > max_diag2 {
                max_diag2 = j1;
                imax_diag = j;
            }
        }
        j += 1;
        k += 1;
    }
    let mlen = if imax_diag > len1 { nall - imax_diag } else { imax_diag };
    let emax = ((1.0 - options.cluster_thd) * mlen as f64) as i32 + 1;
    for j in from..imax_diag {
        // if aap pairs fail to open gap
        if (imax_diag - j) > emax || unpack_diag(ds[j as usize]).0 < 1 {
            best_score -= unpack_diag(ds[j as usize]).0;
            from += 1;
        } else {
            break;
        }
    }
    let mut j = end;
    while j > imax_diag {
        // if aap pairs fail to open gap
        if (j - imax_diag) > emax || unpack_diag(ds[j as usize]).0 < 1 {
            best_score -= unpack_diag(ds[j as usize]).0;
            end -= 1;
        } else {
            break;
        }
        j -= 1;
    }

    *band_left = from - len1 + 1;
    *band_right = end - len1 + 1;
    *band_center = imax_diag - len1 + 1;
    *best_sum = best_score;
    if options.is_454 {
        if *band_left > 0 {
            *best_sum = 0;
        }
        if *band_right < 0 {
            *best_sum = 0;
        }
    }
    OK_FUNC
}
// END diag_test_aapn_est

/*
local alignment of two sequence within a diag band
for band 0 means direction (0,0) -> (1,1)
         1 means direction (0,1) -> (1,2)
        -1 means direction (1,0) -> (2,1)
added on 2006 11 13
band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                            XXXXXXXXXXXXXXX                  seq1
band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                             XXXXXXXXXXXXXXX                 seq1
extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
    band = len2-1                            XXXXXXXXXXXXXXX seq1
band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                           XXXXXXXXXXXXXXX                   seq1
extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
*/
/// `local_band_align`
///
/// The upstream `PRINT` and `MAKEALIGN` blocks are `#ifdef`-guarded on macros that are
/// never defined (each block re-`#define`s its own guard from inside), so they are dead in
/// every shipped build and are not reproduced.
pub fn local_band_align(
    iseq1: &[u8],
    iseq2: &[u8],
    len1: i32,
    len2: i32,
    mat: &ScoreMatrix,
    best_score: &mut i32,
    iden_no: &mut i32,
    alnln: &mut i32,
    dist: &mut f32,
    alninfo: Option<&mut [i32; 5]>,
    band_left: i32,
    band_center: i32,
    band_right: i32,
    buffer: &mut WorkingBuffer,
    options: &Options,
) -> i32 {
    *iden_no = 0;

    if band_right >= len2 || band_left <= -len1 || band_left > band_right {
        return FAILED_FUNC;
    }

    // allocate mem for score_mat[len1][len2] etc
    let band_width = band_right - band_left + 1;
    let band_width1 = (band_width + 1) as usize;

    // score_mat, back_mat [i][j]: i index of seqi (0 to len(seqi)-1), j index of band (0 to band_width-1)
    let rows = len1 as usize + 1;
    // One column of slack past `band_width1`: the sweep below slices `prev_s[j1 + 1]` for
    // the whole `j1` range, and the last `j1` can be `band_width` itself. The extra column
    // is never read -- the `j1 + 1 < band_width` guard rules it out -- it just lets the
    // slice be formed without a special case.
    if band_width1 + 1 > buffer.dp_stride {
        buffer.dp_stride = band_width1 + 1;
        buffer.score_mat.clear();
        buffer.back_mat.clear();
    }
    let stride = buffer.dp_stride;
    let need = rows * stride;
    if buffer.score_mat.len() < need {
        buffer.score_mat.resize(need, 0);
        buffer.back_mat.resize(need, 0);
    }
    let score_mat = &mut buffer.score_mat;
    let back_mat = &mut buffer.back_mat;

    *best_score = 0;
    /* seq1 is query, seq2 is rep
                  seq2    len2 = 17       seq2    len2 = 17    seq2    len2 = 17
                  01234567890123456       01234567890123456    01234567890123456
       0          xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx    xXXXXXXXxxxxxxxxx
       1     \\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx    xx\xxxxx\xxxxxxxx
       2      \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx    xxx\xxxxx\xxxxxxx
  seq1 3       \\\Xx\xxxxxxxxxxxxxx    \\\Xxxxx\xxxxxxxxxxx    xxxx\xxxxx\xxxxxx
  len1 4        \\Xxx\xxxxxxxxxxxxx     \\Xxxxxx\xxxxxxxxxx    xxxxx\xxxxx\xxxxx
  = 11 5         \Xxxx\xxxxxxxxxxxx      \Xxxxxxx\xxxxxxxxx    xxxxxx\xxxxx\xxxx
       6          Xxxxx\xxxxxxxxxxx       Xxxxxxxx\xxxxxxxx    xxxxxxx\xxxxx\xxx
       7          x\xxxx\xxxxxxxxxx       x\xxxxxxx\xxxxxxx    xxxxxxxx\xxxxx\xx
       8          xx\xxxx\xxxxxxxxx       xx\xxxxxxx\xxxxxx    xxxxxxxxx\xxxxx\x
       9          xxx\xxxx\xxxxxxxx       xxx\xxxxxxx\xxxxx    xxxxxxxxxx\xxxxx\
       0          xxxx\xxxx\xxxxxxx       xxxx\xxxxxxx\xxxx    xxxxxxxxxxx\xxxxx
                  band_left < 0 (-6)      band_left < 0 (-6)   band_left >=0 (1)
                  band_right < 0 (-1)     band_right >=0 (2)   band_right >=0(7)
                  band_width 6            band_width 9         band_width 7
       init score_mat, and iden_mat (place with upper 'X')
     */

    if band_left < 0 {
        //set score to left border of the matrix within band
        let tband = if band_right < 0 { band_right } else { 0 };
        //for (k=band_left; k<tband; k++)
        for k in band_left..=tband {
            // fixed on 2006 11 14
            let i = -k;
            let j1 = k - band_left;
            // penalty for leading gap opening = penalty for gap extension
            // each of the left side query hunging residues give ext_gap (-1)
            score_mat[i as usize * stride + j1 as usize] = (mat.ext_gap as i64) * i as i64;
            back_mat[i as usize * stride + j1 as usize] = DP_BACK_TOP;
        }
        back_mat[(-tband) as usize * stride + (tband - band_left) as usize] = DP_BACK_NONE;
    }

    if band_right >= 0 {
        //set score to top border of the matrix within band
        let tband = if band_left > 0 { band_left } else { 0 };
        for j in tband..=band_right {
            let j1 = j - band_left;
            score_mat[j1 as usize] = (mat.ext_gap as i64) * j as i64;
            back_mat[j1 as usize] = DP_BACK_LEFT;
        }
        back_mat[(tband - band_left) as usize] = DP_BACK_NONE;
    }

    let gap_open = [mat.gap, mat.ext_gap];
    let max_diag = band_center - band_left;
    let extra_score = [4i64, 3, 2, 1];
    for i in 1..=len1 {
        let mut j0 = 1 - band_left - i;
        let mut j1max = len2 - band_left - i;
        if j0 < 0 {
            j0 = 0;
        }
        if j1max >= band_width {
            j1max = band_width;
        }
        if j1max < j0 {
            continue;
        }
        let j0u = j0 as usize;
        let n = (j1max - j0 + 1) as usize;

        // `iseq1[i-1]`, and hence its substitution-matrix row, is fixed across the sweep.
        let ci = iseq1[(i - 1) as usize] as usize;
        let mrow = &mat.matrix[ci];

        // Split the two rows the sweep touches out of the flat matrices, then narrow each
        // one to exactly the `j1` range. Indexing those by the loop counter lets the bounds
        // checks fall away: this loop runs once per DP cell and is a quarter of the
        // program's time.
        let cur_base = i as usize * stride;
        let prev_base = cur_base - stride;
        let (s_head, s_tail) = score_mat.split_at_mut(cur_base);
        let (b_head, b_tail) = back_mat.split_at_mut(cur_base);

        // The `j1 - 1` cell is whatever the previous iteration wrote, so carry it in a
        // register; only the first iteration has to read it back, from the band border.
        let mut left_s: i64 = if j0u > 0 { s_tail[j0u - 1] } else { 0 };
        let mut left_b: i32 = if j0u > 0 { b_tail[j0u - 1] } else { 0 };

        let ps_cur = &s_head[prev_base + j0u..][..n];
        let ps_next = &s_head[prev_base + j0u + 1..][..n];
        let pb_next = &b_head[prev_base + j0u + 1..][..n];
        let cs_out = &mut s_tail[j0u..][..n];
        let cb_out = &mut b_tail[j0u..][..n];
        let seq2 = &iseq2[(j0 + i + band_left - 1) as usize..][..n];

        for k in 0..n {
            let j1 = j0 + k as i32;
            let j = j1 + i + band_left;

            let cj = seq2[k] as usize;
            let mut sij = mrow[cj] as i64;
            //int iden_ij = (ci == cj);

            /* extra score according to the distance to the best diagonal */
            let extra = extra_score[((j1 - max_diag).abs() & 3) as usize]; // max distance 3
            sij += extra * (sij > 0) as i64;

            let mut back = DP_BACK_LEFT_TOP;
            let mut best_score1 = ps_cur[k] + sij;
            let gap0 = gap_open[((i == len1) | (j == len2)) as usize];
            let mut gap;

            if j1 > 0 {
                gap = gap0;
                if left_b == DP_BACK_LEFT {
                    gap = mat.ext_gap;
                }
                let score = left_s + gap as i64;
                if score > best_score1 {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            if j1 + 1 < band_width {
                gap = gap0;
                if pb_next[k] == DP_BACK_TOP {
                    gap = mat.ext_gap;
                }
                let score = ps_next[k] + gap as i64;
                if score > best_score1 {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            cs_out[k] = best_score1;
            cb_out[k] = back;
            left_s = best_score1;
            left_b = back;
        }
    }
    let mut i;
    let mut j;
    if len2 - band_left < len1 {
        i = len2 - band_left;
        j = len2;
    } else if len1 + band_right < len2 {
        i = len1;
        j = len1 + band_right;
    } else {
        i = len1;
        j = len2;
    }
    let mut j1 = j - i - band_left;
    // COMPAT: upstream narrows the int64 cell into the `int&` out-parameter. Scores are
    // scaled by MAX_SEQ, so this genuinely truncates on long alignments; preserved.
    *best_score = score_mat[i as usize * stride + j1 as usize] as i32;
    let best_score1 = score_mat[i as usize * stride + j1 as usize];

    let has_alninfo = alninfo.is_some();
    let mut back = back_mat[i as usize * stride + j1 as usize];
    let mut last = back;
    let mut count = 0i32;
    let mut count2 = 0i32;
    let mut count3 = 0i32;
    let mut begin1 = 0i32;
    let mut begin2 = 0i32;
    let mut end1 = 0i32;
    let mut end2 = 0i32;
    let mut gbegin1 = 0i32;
    let mut gbegin2 = 0i32;
    let mut gend1 = 0i32;
    let mut gend2 = 0i32;
    let mut smin = best_score1;
    let mut smax = best_score1 - 1;
    let mut posmin = 0i32;
    let mut posmax = 0i32;
    let mut pos = 0i32;
    let mut dlen = 0i32;
    let mut dcount = 0i32;

    let mut masked = 0i32;
    let mut indels = 0i32;
    let mut max_indels = 0i32;
    while back != DP_BACK_NONE {
        match back {
            DP_BACK_TOP => {
                let bl = ((last != back) & (j != 1) & (j != len2)) as i32;
                dlen += bl;
                dcount += bl;
                let score = score_mat[i as usize * stride + j1 as usize];
                if score < smin {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                j1 += 1;
            }
            DP_BACK_LEFT => {
                let bl = ((last != back) & (i != 1) & (i != len1)) as i32;
                dlen += bl;
                dcount += bl;
                let score = score_mat[i as usize * stride + j1 as usize];
                if score < smin {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j1 -= 1;
                j -= 1;
            }
            DP_BACK_LEFT_TOP => {
                if has_alninfo && options.global_identity {
                    if i == 1 || j == 1 {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if i == len1 || j == len2 {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                let score = score_mat[i as usize * stride + j1 as usize];
                i -= 1;
                j -= 1;
                let is_match = iseq1[i as usize] == iseq2[j as usize];
                if score > smax {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if options.is_est && (iseq1[i as usize] > 4 || iseq2[j as usize] > 4) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !is_match as i32;
                    count += is_match as i32;
                    count2 += is_match as i32;
                    count3 += is_match as i32;
                }
                if score < smin {
                    let mm = !is_match as i32;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
            }
            _ => println!("{}", back),
        }
        if options.is_454 {
            if back == DP_BACK_LEFT_TOP {
                if indels > max_indels {
                    max_indels = indels;
                }
                indels = 0;
            } else if last == DP_BACK_LEFT_TOP {
                indels = 1;
            } else if indels != 0 {
                indels += 1;
            }
        }
        pos += 1;
        last = back;
        back = back_mat[i as usize * stride + j1 as usize];
    }
    if options.is_454 && max_indels > options.max_indel {
        return FAILED_FUNC;
    }
    *iden_no = if options.global_identity {
        count3
    } else {
        count - count2
    };
    *alnln = posmin - posmax + 1 - masked;
    *dist = dcount as f32 / dlen as f32;
    //dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
    let umtail1 = len1 - 1 - end1;
    let umtail2 = len2 - 1 - end2;
    let umhead = if begin1 < begin2 { begin1 } else { begin2 };
    let umtail = if umtail1 < umtail2 { umtail1 } else { umtail2 };
    let umlen = umhead + umtail;
    if umlen > options.unmatch_len {
        return FAILED_FUNC;
    }
    if (umlen as f64) > len1 as f64 * options.short_unmatch_per {
        return FAILED_FUNC;
    }
    if (umlen as f64) > len2 as f64 * options.long_unmatch_per {
        return FAILED_FUNC;
    }
    if let Some(info) = alninfo {
        info[0] = begin1;
        info[1] = end1;
        info[2] = begin2;
        info[3] = end2;
        info[4] = masked;
        if options.global_identity {
            info[0] = gbegin1;
            info[1] = gend1;
            info[2] = gbegin2;
            info[3] = gend2;
        }
    }

    OK_FUNC
}
// END int local_band_align

impl Alphabet {
    /// The default protein alphabet, i.e. the initial values of `aa2idx` and `na2idx`.
    pub fn new() -> Alphabet {
        Alphabet {
            aa2idx: AA2IDX,
            na2idx: NA2IDX,
        }
    }

    /// `setaa_to_na`
    pub fn setaa_to_na(&mut self) {
        for i in 0..26 {
            self.aa2idx[i] = self.na2idx[i];
        }
    }
}

impl Default for Alphabet {
    fn default() -> Self {
        Alphabet::new()
    }
}

/////////////////
impl ScoreMatrix {
    /// `ScoreMatrix::ScoreMatrix`
    pub fn new() -> ScoreMatrix {
        let mut m = ScoreMatrix {
            matrix: [[0; MAX_AA]; MAX_AA],
            gap: 0,
            ext_gap: 0,
        };
        m.init();
        m
    }

    /// `ScoreMatrix::init`
    pub fn init(&mut self) {
        self.set_gap(-11, -1);
        self.set_matrix(&BLOSUM62);
    }

    /// `ScoreMatrix::set_gap`
    pub fn set_gap(&mut self, gap1: i32, ext_gap1: i32) {
        self.gap = MAX_SEQ.wrapping_mul(gap1);
        self.ext_gap = MAX_SEQ.wrapping_mul(ext_gap1);
    }

    /// `ScoreMatrix::set_matrix`: read `mat1` as a lower triangle and mirror it.
    ///
    /// COMPAT: upstream always reads `MAX_AA * (MAX_AA + 1) / 2 == 276` entries, which
    /// overruns `BLOSUM62_na` (28 entries) by 248 ints. Only cells with both indices
    /// `<= 5` are ever looked up for nucleotides, and those come from `k <= 20`, so the
    /// overread never reaches a cell that is used; the out-of-range reads yield 0 here
    /// instead of whatever followed the array in the C++ binary's `.data`.
    pub fn set_matrix(&mut self, mat1: &[i32]) {
        let mut k = 0;
        for i in 0..MAX_AA {
            for j in 0..=i {
                let v = MAX_SEQ.wrapping_mul(mat1.get(k).copied().unwrap_or(0));
                k += 1;
                self.matrix[j][i] = v;
                self.matrix[i][j] = v;
            }
        }
    }

    /// `ScoreMatrix::set_to_na`
    pub fn set_to_na(&mut self) {
        self.set_gap(-6, -1);
        self.set_matrix(&BLOSUM62_NA);
    }

    /// `ScoreMatrix::set_match` -- only for est.
    pub fn set_match(&mut self, score: i32) {
        for i in 0..5 {
            self.matrix[i][i] = MAX_SEQ.wrapping_mul(score);
        }
        //matrix[3][4] = matrix[4][3] = MAX_SEQ * score;
    }

    /// `ScoreMatrix::set_mismatch` -- only for est.
    pub fn set_mismatch(&mut self, score: i32) {
        for i in 0..MAX_AA {
            for j in 0..i {
                let v = MAX_SEQ.wrapping_mul(score);
                self.matrix[j][i] = v;
                self.matrix[i][j] = v;
            }
        }
        self.matrix[3][4] = MAX_SEQ;
        self.matrix[4][3] = MAX_SEQ;
    }
}

impl Default for ScoreMatrix {
    fn default() -> Self {
        ScoreMatrix::new()
    }
}

impl WordTable {
    /// `WordTable::WordTable`
    pub fn new(naa: i32, naan: i32) -> WordTable {
        let mut t = WordTable {
            row_head: Vec::new(),
            arena: Vec::new(),
            free: Vec::new(),
            sequences: Vec::new(),
            naa: 0,
            naan: 0,
            is_aa: 1,
            size: 0,
            frag_count: 0,
        };
        t.init(naa, naan);
        t
    }

    /// `WordTable::SetDNA`
    ///
    /// COMPAT: never called upstream; `is_aa` is therefore always 1 and is never read.
    pub fn set_dna(&mut self) {
        self.is_aa = 0;
    }

    /// `WordTable::Init`
    pub fn init(&mut self, naa: i32, naan: i32) {
        self.naa = naa;
        self.naan = naan;
        self.row_head.clear();
        self.row_head.resize(naan as usize, NO_BLOCK);
    }

    /// `WordTable::Clear`
    ///
    /// Upstream sets each row's logical size to 0 without freeing, so the rows keep their
    /// capacity for the next block; `Vec::clear` has the same effect.
    pub fn clear(&mut self) {
        self.size = 0;
        self.frag_count = 0;
        self.sequences.clear();
        self.row_head.fill(NO_BLOCK);
        self.arena.clear();
        for list in self.free.iter_mut() {
            list.clear();
        }
    }

    /// Append one entry to word `j`, mirroring `indexCounts[j].Append(ic)`.
    #[inline]
    fn push_entry(&mut self, j: usize, ic: IndexCount) {
        let off = self.row_head[j];
        if off == NO_BLOCK {
            let o = self.alloc_block(0) as usize;
            self.arena[o].count = 1;
            self.arena[o + 1] = ic;
            self.row_head[j] = o as u32;
            return;
        }
        let o = off as usize;
        let cap = self.arena[o].index as u32;
        let len = self.arena[o].count as u32;
        if len < cap {
            self.arena[o + 1 + len as usize] = ic;
            self.arena[o].count = (len + 1) as i32;
            return;
        }
        // Full: move the row to the next size class and hand the old block back.
        let k = word_block_class(cap);
        let new_off = self.alloc_block(k + 1) as usize;
        self.arena
            .copy_within(o + 1..o + 1 + len as usize, new_off + 1);
        self.arena[new_off + 1 + len as usize] = ic;
        self.arena[new_off].count = (len + 1) as i32;
        self.free[k as usize].push(off);
        self.row_head[j] = new_off as u32;
    }

    /// Take a block of size class `k`, from the free list if one is waiting.
    #[inline]
    fn alloc_block(&mut self, k: u32) -> u32 {
        if self.free.len() <= k as usize {
            self.free.resize(k as usize + 1, Vec::new());
        }
        if let Some(off) = self.free[k as usize].pop() {
            self.arena[off as usize].count = 0;
            return off;
        }
        let cap = word_block_cap(k);
        if self.arena.len() + 1 + cap as usize >= NO_BLOCK as usize {
            // Offsets are 32-bit; this needs a word table well past 30 GB, which `-M`
            // normally rules out, but say so rather than wrapping.
            bomb_error("word table too large, please lower -M");
        }
        let off = self.arena.len() as u32;
        self.arena.push(IndexCount::new(cap as i32, 0));
        self.arena
            .resize(off as usize + 1 + cap as usize, IndexCount::default());
        off
    }

    // COMPAT: the `AddWordCounts(NVector<IndexCount>&, Sequence*, bool)` and
    // `AddWordCountsFrag(NVector<IndexCount>&, int, int, int)` overloads have no callers
    // upstream (the latter is a stub that returns 0). Not reproduced.

    /// `WordTable::AddWordCounts(int, Vector<int>&, Vector<INTs>&, int, bool)`
    pub fn add_word_counts(
        &mut self,
        aan_no: i32,
        word_encodes: &[i32],
        word_encodes_no: &[INTs],
        idx: i32,
        skip_n: bool,
    ) -> i32 {
        for i in 0..aan_no as usize {
            let k = word_encodes_no[i];
            if k != 0 {
                let j = word_encodes[i];
                if skip_n && j < 0 {
                    continue; // for those has 'N'
                }
                self.push_entry(j as usize, IndexCount::new(idx, k as i32));
                self.size += 1;
            }
        }
        OK_FUNC
    }

    /// `WordTable::AddWordCountsFrag(int, Vector<int>&, Vector<INTs>&, int, int)`
    ///
    /// Takes `word_encodes` by value-like mutable reference because it sorts each fragment
    /// of the array in place, exactly as upstream does.
    pub fn add_word_counts_frag(
        &mut self,
        aan_no: i32,
        word_encodes: &mut [i32],
        word_encodes_no: &mut [INTs],
        frag: i32,
        frag_size: i32,
    ) -> i32 {
        for i in 0..frag {
            let k = if (i + 1) * frag_size < aan_no {
                (i + 1) * frag_size - 1
            } else {
                aan_no - 1
            };
            //quick_sort(&word_encodes[0], i*frag_size, k);
            word_encodes[(i * frag_size) as usize..(k + 1) as usize].sort_unstable();
        }
        let mut j = aan_no - 1;
        while j > 0 {
            if word_encodes[j as usize] == word_encodes[(j - 1) as usize] {
                word_encodes_no[(j - 1) as usize] += word_encodes_no[j as usize];
                word_encodes_no[j as usize] = 0;
            }
            j -= 1;
        }
        // END check_word_encodes

        let mut i = 0;
        while i < aan_no {
            let k = if frag_size < (aan_no - i) {
                frag_size
            } else {
                aan_no - i
            };
            let fra = i / frag_size;
            //AddWordCounts(k, word_encodes+i, word_encodes_no+i, NR90f_no+fra);
            for i1 in i..i + k {
                let k1 = word_encodes_no[i1 as usize];
                if k1 != 0 {
                    let j = word_encodes[i1 as usize];
                    self.push_entry(j as usize, IndexCount::new(self.frag_count + fra, k1 as i32));
                    self.size += 1;
                }
            }
            i += frag_size;
        }
        self.frag_count += frag;

        0
    }

    /// Walk the entries of word `j` in insertion order.
    fn row_iter(&self, j: usize) -> impl Iterator<Item = IndexCount> + '_ {
        let off = self.row_head[j];
        let (start, len) = if off == NO_BLOCK {
            (0usize, 0usize)
        } else {
            (off as usize + 1, self.arena[off as usize].count as usize)
        };
        self.arena[start..start + len].iter().copied()
    }

    /// `WordTable::PrintAll`
    pub fn print_all(&self) {
        let mut cols = 0;
        let mut total_words: i64 = 0;
        let mut k = 0;
        for i in 0..self.naan as usize {
            let size: usize = self.row_iter(i).count();
            if size == 0 {
                continue;
            }
            cols += 1;
            print!("{}\t{}\tsize:{}\t", k, i, size);
            for ic in self.row_iter(i) {
                print!("{},{} ", ic.index, ic.count);
                total_words += ic.count as i64;
            }
            println!();
            k += 1;
        }

        println!("total cols: {} total words: {}", cols, total_words);
    }

    /// `WordTable::CountWords`
    ///
    /// Accumulates, for every representative in the table, how many words it shares with
    /// the query. `look_counts` collects one entry per representative that was hit;
    /// `index_mapping[rep]` is that entry's 1-based position, or 0 when unseen.
    pub fn count_words(
        &self,
        aan_no: i32,
        word_encodes: &[i32],
        word_encodes_no: &[INTs],
        look_counts: &mut LookCounts,
        index_mapping: &mut [u32],
        est: bool,
        min: i32,
    ) -> i32 {
        for j in 0..look_counts.size {
            index_mapping[look_counts.items[j].index as usize] = 0;
        }
        look_counts.size = 0;

        // Hoist the two hot buffers out of their owners: the inner loop runs once per
        // (query word, representative) pair and is the single hottest loop in the program,
        // so re-loading `look_counts.items`' pointer through the struct on every hit is
        // worth avoiding.
        let items = &mut look_counts.items[..];
        let mut size = look_counts.size;

        let we = word_encodes;
        let wen = word_encodes_no;
        let mut j0 = 0usize;
        if est {
            // COMPAT: upstream scans `while( *we < 0 ) j0++, we++` with no bound. For an
            // all-N sequence every word in [0, aan_no) is -1, so it runs off into stale
            // entries left by an earlier call; wherever it stops is at or past `aan_no`,
            // and the loop below then does nothing. Bounding the scan by the buffer length
            // reaches the same state without the overread.
            while j0 < we.len() && we[j0] < 0 {
                j0 += 1; // if met short word has 'N'
            }
        }
        for (k, (&j, &n)) in we
            .iter()
            .zip(wen.iter())
            .enumerate()
            .take(aan_no as usize)
            .skip(j0)
        {
            let j1 = n as i32;
            if j1 == 0 {
                continue;
            }
            let rest = aan_no - k as i32 + 1;
            let off = self.row_head[j as usize];
            if off == NO_BLOCK {
                continue;
            }
            let base = off as usize;
            let row_len = self.arena[base].count as usize;
            for ic in &self.arena[base + 1..base + 1 + row_len] {
                let c = if ic.count < j1 { ic.count } else { j1 };
                debug_assert!((ic.index as usize) < index_mapping.len());
                debug_assert!(size + 1 < items.len());
                // SAFETY: `ic.index` is a representative's position in `table.sequences`,
                // which never exceeds `max_sequences`; `index_mapping` and `look_counts`
                // are both sized `frag + 2` for that same bound in `WorkingBuffer::set`.
                // `size` counts distinct representatives seen so far, so it too stays below
                // that bound, and `*idm - 1` is an index this loop wrote itself. This is the
                // innermost loop of the program -- roughly half its run time -- and the
                // three bounds checks are a measurable share of it.
                let idm = unsafe { index_mapping.get_unchecked_mut(ic.index as usize) };
                if *idm == 0 {
                    if rest < min {
                        continue;
                    }
                    unsafe { *items.get_unchecked_mut(size) = IndexCount::new(ic.index, c) };
                    size += 1;
                    *idm = size as u32;
                } else {
                    unsafe { items.get_unchecked_mut((*idm - 1) as usize).count += c };
                }
            }
        }
        look_counts.size = size;
        look_counts.items[look_counts.size].count = 0;
        OK_FUNC
    }
}

/* Quick Sort.
 * Adam Drozdek: Data Structures and Algorithms in C++, 2nd Edition.
 */
/// `PartialQuickSort`
///
/// COMPAT: has no callers upstream; translated for completeness.
pub fn partial_quick_sort(data: &mut [IndexCount], first: i32, last: i32, partial: i32) {
    let mut lower = first + 1;
    let mut upper = last;
    if first >= last {
        return;
    }
    data.swap(first as usize, ((first + last) / 2) as usize);
    let pivot = data[first as usize];

    while lower <= upper {
        while lower <= last && data[lower as usize].count < pivot.count {
            lower += 1;
        }
        while pivot.count < data[upper as usize].count {
            upper -= 1;
        }
        if lower < upper {
            data.swap(lower as usize, upper as usize);
            upper -= 1;
        }
        lower += 1;
    }
    data.swap(first as usize, upper as usize);
    if first < upper - 1 {
        partial_quick_sort(data, first, upper - 1, partial);
    }
    if upper >= partial {
        return;
    }
    if upper + 1 < last {
        partial_quick_sort(data, upper + 1, last, partial);
    }
}

impl LookCounts {
    /// Sized as upstream sizes `lookCounts`: one slot per table entry plus two, so the
    /// terminating zero-count entry always has somewhere to go.
    pub fn new(n: usize) -> LookCounts {
        LookCounts {
            items: vec![IndexCount::default(); n],
            size: 0,
        }
    }
}

impl Default for Sequence {
    /// `Sequence::Sequence()`
    fn default() -> Sequence {
        Sequence {
            data: Vec::new(),
            size: 0,
            size_r2: 0,
            swap: None,
            offset: 0,
            des_begin: 0,
            des_begin2: 0,
            tot_length: 0,
            tot_length2: 0,
            identifier: Vec::new(),
            index: 0,
            state: 0,
            cluster_id: 0,
            identity: 0.0,
            distance: 2.0,
            coverage: [0; 4],
        }
    }
}

impl Clone for Sequence {
    /// `Sequence::Sequence( const Sequence & other )`
    fn clone(&self) -> Sequence {
        Sequence {
            data: self.data.clone(),
            size: self.size,
            size_r2: if self.data.is_empty() { self.size_r2 } else { 0 },
            swap: self.swap,
            offset: self.offset,
            des_begin: self.des_begin,
            des_begin2: self.des_begin2,
            tot_length: self.tot_length,
            tot_length2: self.tot_length2,
            identifier: self.identifier.clone(),
            index: self.index,
            state: self.state,
            cluster_id: self.cluster_id,
            identity: self.identity,
            distance: 2.0,
            coverage: self.coverage,
        }
    }
}

impl Sequence {
    // back to back merge for PE
    // R1 -> XXXXXXABC ------------------- NMLYYYYYY <--R2
    // >R1           >R2
    // XXXXXXABC     YYYYYYLMN =====> Merge into
    // >R12
    // NMLYYYYYYXXXXXXABC
    /// `Sequence::Sequence( const Sequence & other, const Sequence & other2, int mode )`
    pub fn from_pair(other: &Sequence, other2: &Sequence, mode: i32) -> Sequence {
        if mode != 1 {
            bomb_error("unknown mode");
        }

        // Upstream memcpy's all of `other`'s fields, then overrides the ones below.
        let mut s = Sequence {
            data: Vec::new(),
            size: other.size,
            size_r2: other.size_r2,
            swap: other.swap,
            offset: other.offset,
            des_begin: other.des_begin,
            des_begin2: other.des_begin2,
            tot_length: other.tot_length,
            tot_length2: other.tot_length2,
            identifier: other.identifier.clone(), // only use R1
            index: other.index,
            state: other.state,
            cluster_id: other.cluster_id,
            identity: other.identity,
            distance: 2.0,
            coverage: other.coverage,
        };

        if !other.data.is_empty() && !other2.data.is_empty() {
            s.size = other.size + other2.size;
            s.size_r2 = other2.size;
            s.data.extend_from_slice(&other2.data[..s.size_r2 as usize]); // copy R2 first
            strrev(&mut s.data); // reverse R2 on data
            s.data
                .extend_from_slice(&other.data[..(s.size - s.size_r2) as usize]); // copy R1 to end of R2
            s.des_begin2 = other2.des_begin;
            s.tot_length2 = other2.tot_length;
        } else if !other.data.is_empty() || !other2.data.is_empty() {
            bomb_error("Not both PE sequences have data");
        }

        s
    }

    /// Read `data[i]`, yielding the C buffer's NUL terminator for `i == size`.
    ///
    /// The C++ `Sequence` keeps `data` NUL-terminated, and the header parser relies on
    /// being able to read one byte past the last one when the line is empty.
    #[inline]
    pub fn data_at(&self, i: usize) -> u8 {
        self.data.get(i).copied().unwrap_or(0)
    }

    /// `Sequence::Clear`
    ///
    /// Deliberately leaves `size` alone: it is still needed when writing the output.
    pub fn clear(&mut self) {
        self.data = Vec::new();
    }

    /// `Sequence::operator=( const char *s )`
    pub fn assign(&mut self, s: &[u8]) {
        self.data.clear();
        self.data.extend_from_slice(s);
        self.size = self.data.len() as i32;
    }

    /// `Sequence::operator+=( const char *s )`
    pub fn append(&mut self, s: &[u8]) {
        self.data.extend_from_slice(s);
        self.size = self.data.len() as i32;
    }

    /// `Sequence::Resize`
    pub fn resize(&mut self, n: i32) {
        self.data.resize(n as usize, 0);
        self.size = n;
    }

    /// `Sequence::trim`
    pub fn trim(&mut self, trim_len: i32) {
        if trim_len >= self.size {
            return;
        }
        self.size = trim_len;
        self.data.truncate(trim_len as usize);
    }

    /// `Sequence::ConvertBases`
    pub fn convert_bases(&mut self, alphabet: &Alphabet) {
        for i in 0..self.size as usize {
            self.data[i] = alphabet.aa2idx[(self.data[i] - b'A') as usize] as u8;
        }
    }

    /// `Sequence::Format`: returns 0 when the record is usable, otherwise the number of
    /// offending characters.
    pub fn format(&mut self) -> i32 {
        while self.size != 0 && c_isspace(self.data[(self.size - 1) as usize]) {
            self.size -= 1;
        }
        if self.size != 0 && self.data[(self.size - 1) as usize] == b'*' {
            self.size -= 1;
        }
        self.data.truncate(self.size as usize);
        let mut m = 0;
        for i in 0..self.size as usize {
            let ch = self.data[i];
            m += !(ch.is_ascii_alphabetic() || c_isspace(ch)) as i32;
        }
        if m != 0 {
            return m;
        }
        let mut j = 0usize;
        for i in 0..self.size as usize {
            let ch = self.data[i];
            if ch.is_ascii_alphabetic() {
                self.data[j] = ch.to_ascii_uppercase();
                j += 1;
            }
        }
        self.data.truncate(j);
        self.size = j as i32;
        0
    }

    // COMPAT: `SwapIn` / `SwapOut` are inert upstream -- nothing ever assigns
    // `Sequence::swap`, so `-B` has no effect. Not reproduced; see the TempFile note.

    /// `Sequence::PrintInfo`
    pub fn print_info<W: std::io::Write>(&self, id: i32, fout: &mut W, options: &Options) {
        let tag = if options.is_est { "nt" } else { "aa" };
        let print = options.print != 0;
        let strand = options.is_est;
        let _ = write!(
            fout,
            "{}\t{}{}, >{}...",
            id,
            self.size,
            tag,
            String::from_utf8_lossy(&self.identifier[1..])
        );
        if self.identity != 0.0 {
            let c = &self.coverage;
            let _ = write!(fout, " at ");
            if print {
                let _ = write!(fout, "{}:{}:{}:{}/", c[0], c[1], c[2], c[3]);
            }
            if strand {
                let _ = write!(
                    fout,
                    "{}/",
                    if self.state & IS_MINUS_STRAND != 0 {
                        '-'
                    } else {
                        '+'
                    }
                );
            }
            let _ = write!(fout, "{:.2}%", self.identity * 100.0);
            if options.use_distance {
                let _ = write!(fout, "/{:.2}%", self.distance * 100.0);
            }
            let _ = writeln!(fout);
        } else {
            let _ = writeln!(fout, " *");
        }
    }
}

/// Byte-oriented `fgets` over a reader, tracking the (uncompressed) stream offset so that
/// `ftell` / `gztell` can be reproduced for `Sequence::des_begin`.
pub struct LineReader<R: std::io::Read> {
    inner: std::io::BufReader<R>,
    pos: u64,
}

impl<R: std::io::Read> LineReader<R> {
    pub fn new(inner: R) -> LineReader<R> {
        LineReader {
            inner: std::io::BufReader::with_capacity(1 << 16, inner),
            pos: 0,
        }
    }

    /// `fgets( buf, max, stream )`: read up to `max - 1` bytes, stopping after a newline.
    /// Returns false at end of input with nothing read, in which case `buf` is untouched --
    /// exactly like C, which leaves the caller's buffer alone when it returns NULL.
    pub fn fgets(&mut self, buf: &mut Vec<u8>, max: usize) -> bool {
        use std::io::BufRead;
        let limit = max - 1;
        let mut got = 0usize;
        let mut out: Vec<u8> = Vec::new();
        loop {
            let available = match self.inner.fill_buf() {
                Ok(b) => b,
                Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(_) => break,
            };
            if available.is_empty() {
                break;
            }
            let room = limit - got;
            let take = available.len().min(room);
            match available[..take].iter().position(|&c| c == b'\n') {
                Some(p) => {
                    out.extend_from_slice(&available[..=p]);
                    self.inner.consume(p + 1);
                    self.pos += (p + 1) as u64;
                    got += p + 1;
                    break;
                }
                None => {
                    out.extend_from_slice(&available[..take]);
                    self.inner.consume(take);
                    self.pos += take as u64;
                    got += take;
                    if got >= limit {
                        break;
                    }
                }
            }
        }
        if got == 0 {
            return false;
        }
        buf.clear();
        buf.extend_from_slice(&out);
        true
    }

    /// `ftell` / `gztell`
    pub fn tell(&self) -> u64 {
        self.pos
    }
}

/// The input database reopened for output writing: `WriteClusters` seeks back to each
/// representative's `des_begin` and copies `tot_length` bytes verbatim.
///
/// Plain files seek directly. Gzip streams cannot seek backwards, so a backward seek
/// restarts decompression from the beginning -- which is what `gzseek` does internally,
/// and why the `-sf` help text warns that it is slow on `.gz` input.
enum RecordSource {
    Plain(std::fs::File),
    Gz {
        path: String,
        inner: flate2::read::MultiGzDecoder<std::fs::File>,
        pos: u64,
    },
}

impl RecordSource {
    fn open(path: &str) -> RecordSource {
        let f = match std::fs::File::open(path) {
            Ok(f) => f,
            Err(_) => bomb_error("file opening failed"),
        };
        if is_gz(path) {
            RecordSource::Gz {
                path: path.to_string(),
                inner: flate2::read::MultiGzDecoder::new(f),
                pos: 0,
            }
        } else {
            RecordSource::Plain(f)
        }
    }

    /// `fseek( fin, offset, SEEK_SET )` / `gzseek( fin, offset, SEEK_SET )`
    fn seek(&mut self, offset: u64) {
        match self {
            RecordSource::Plain(f) => {
                let _ = std::io::Seek::seek(f, std::io::SeekFrom::Start(offset));
            }
            RecordSource::Gz { path, inner, pos } => {
                if offset < *pos {
                    let f = match std::fs::File::open(&*path) {
                        Ok(f) => f,
                        Err(_) => bomb_error("file opening failed"),
                    };
                    *inner = flate2::read::MultiGzDecoder::new(f);
                    *pos = 0;
                }
                let mut skip = [0u8; 1 << 16];
                while *pos < offset {
                    let want = ((offset - *pos) as usize).min(skip.len());
                    match std::io::Read::read(inner, &mut skip[..want]) {
                        Ok(0) | Err(_) => break,
                        Ok(n) => *pos += n as u64,
                    }
                }
            }
        }
    }
}

impl std::io::Read for RecordSource {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            RecordSource::Plain(f) => f.read(buf),
            RecordSource::Gz { inner, pos, .. } => {
                let n = inner.read(buf)?;
                *pos += n as u64;
                Ok(n)
            }
        }
    }
}

/// The body of `WriteClusters`: for each selected record, seek to `des_begin` and copy
/// `tot_length` bytes, in `MAX_LINE_SIZE` chunks exactly as upstream does.
fn copy_records<W: std::io::Write>(
    src: &mut RecordSource,
    out: &mut W,
    records: &[(u64, i32)],
) {
    let mut buf = vec![0u8; MAX_LINE_SIZE];
    for &(des_begin, tot_length) in records {
        src.seek(des_begin);
        let count = tot_length as usize / MAX_LINE_SIZE;
        let rest = tot_length as usize % MAX_LINE_SIZE;
        for _ in 0..count {
            if read_exactly(src, &mut buf[..MAX_LINE_SIZE]) == 0 {
                bomb_error("Can not swap in sequence");
            }
            let _ = out.write_all(&buf[..MAX_LINE_SIZE]);
        }
        if rest != 0 {
            if read_exactly(src, &mut buf[..rest]) == 0 {
                bomb_error("Can not swap in sequence");
            }
            let _ = out.write_all(&buf[..rest]);
        }
    }
}

/// `fread`: fill the whole buffer if possible, returning how many bytes were read.
fn read_exactly<R: std::io::Read>(src: &mut R, buf: &mut [u8]) -> usize {
    let mut got = 0;
    while got < buf.len() {
        match src.read(&mut buf[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(_) => break,
        }
    }
    got
}

/// `strcmp(file + strlen(file) - 3, ".gz") == 0`, on bytes so a non-ASCII filename cannot
/// panic on a character boundary.
fn is_gz(path: &str) -> bool {
    path.as_bytes().ends_with(b".gz")
}

fn open_read(path: &str) -> Box<dyn std::io::Read> {
    let f = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(_) => bomb_error("Failed to open the database file"),
    };
    if is_gz(path) {
        Box::new(flate2::read::MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    }
}

impl SequenceDB {
    /// `SequenceDB::Clear`
    pub fn clear(&mut self) {
        self.sequences.clear();
        self.rep_seqs.clear();
    }

    /// `SequenceDB::Read` / `SequenceDB::Readgz`
    ///
    /// The two upstream functions are identical apart from `fgets`/`ftell` versus
    /// `gzgets`/`gztell`, so they share one implementation here parameterised by the
    /// reader; `.gz` is still selected by suffix exactly as upstream does.
    pub fn read(&mut self, file: &str, options: &Options) {
        let mut src = LineReader::new(open_read(file));
        self.read_impl(&mut src, file, options);
    }

    fn read_impl<R: std::io::Read>(
        &mut self,
        fin: &mut LineReader<R>,
        file: &str,
        options: &Options,
    ) {
        let mut one = Sequence::default();
        let mut des = Sequence::default();
        let option_l = options.min_length;
        self.clear();
        let mut buffer: Vec<u8> = vec![b'>'];
        let mut res;

        loop {
            /* do not break when the last sequence is not handled */
            if buffer.is_empty() {
                buffer.push(b'>');
            } else {
                buffer[0] = b'>';
            }
            res = fin.fgets(&mut buffer, MAX_LINE_SIZE);
            if !res && one.size == 0 {
                break;
            }
            if buffer[0] == b'+' {
                let mut len = buffer.len();
                let mut len2 = len;
                while len2 != 0 && buffer[len2 - 1] != b'\n' {
                    if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                        break;
                    }
                    len2 = buffer.len();
                    len += len2;
                }
                one.tot_length += len as i32;

                // read next line quality score
                if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                    bomb_error("can not read quality score after");
                }
                len = buffer.len();
                len2 = len;
                while len2 != 0 && buffer[len2 - 1] != b'\n' {
                    if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                        break;
                    }
                    len2 = buffer.len();
                    len += len2;
                }
                one.tot_length += len as i32;
            } else if buffer[0] == b'>' || buffer[0] == b'@' || (!res && one.size != 0) {
                if one.size != 0 {
                    // write previous record
                    if one.identifier.is_empty() || one.format() != 0 {
                        println!("Warning: from file \"{}\",", file);
                        println!(
                            "Discarding invalid sequence or sequence without identifier and description!\n"
                        );
                        if !one.identifier.is_empty() {
                            println!("{}", String::from_utf8_lossy(&one.identifier));
                        }
                        println!("{}", String::from_utf8_lossy(&one.data));
                        one.size = 0;
                        one.data.clear();
                    }
                    one.index = self.sequences.len() as i32;
                    if one.size > option_l {
                        if options.trim_len > 0 {
                            one.trim(options.trim_len);
                        }
                        self.sequences.push(one.clone());
                    }
                }
                one.size = 0;
                one.data.clear();
                one.tot_length = 0;

                let mut len = buffer.len();
                let mut len2 = len;
                des.size = 0;
                des.data.clear();
                des.append(&buffer);
                while len2 != 0 && buffer[len2 - 1] != b'\n' {
                    if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                        break;
                    }
                    des.append(&buffer);
                    len2 = buffer.len();
                    len += len2;
                }
                let offset = fin.tell();
                one.des_begin = offset - len as u64;
                one.tot_length += len as i32; // count first line

                let mut i = 0usize;
                if des.data_at(i) == b'>' || des.data_at(i) == b'@' || des.data_at(i) == b'+' {
                    i += 1;
                }
                if des.data_at(i) == b' ' || des.data_at(i) == b'\t' {
                    i += 1;
                }
                if options.des_len != 0 && options.des_len < des.size {
                    des.size = options.des_len;
                }
                while i < des.size as usize && !c_isspace(des.data_at(i)) {
                    i += 1;
                }
                one.identifier.clear();
                one.identifier
                    .extend_from_slice(&des.data[..i.min(des.data.len())]);
            } else {
                one.tot_length += buffer.len() as i32;
                one.append(&buffer);
            }
        }
    }

    /// `SequenceDB::Read( db, db_pe, options )` / `SequenceDB::Readgz( db, db_pe, options )`
    pub fn read_pe(&mut self, file: &str, file2: &str, options: &Options) {
        if is_gz(file) && !is_gz(file2) {
            bomb_error("Both input files need to be in .gz format");
        }
        let mut fin = LineReader::new(open_read(file));
        let mut fin2 = LineReader::new(open_read(file2));
        self.read_pe_impl(&mut fin, &mut fin2, file, file2, options);
    }

    fn read_pe_impl<R: std::io::Read>(
        &mut self,
        fin: &mut LineReader<R>,
        fin2: &mut LineReader<R>,
        file: &str,
        file2: &str,
        options: &Options,
    ) {
        let mut one = Sequence::default();
        let mut two = Sequence::default();
        let mut des = Sequence::default();
        let option_l = options.min_length;
        self.clear();
        let mut buffer: Vec<u8> = vec![b'>'];
        let mut buffer2: Vec<u8> = vec![b'>'];

        loop {
            /* do not break when the last sequence is not handled */
            if buffer.is_empty() {
                buffer.push(b'>');
            } else {
                buffer[0] = b'>';
            }
            let res = fin.fgets(&mut buffer, MAX_LINE_SIZE);
            if buffer2.is_empty() {
                buffer2.push(b'>');
            } else {
                buffer2[0] = b'>';
            }
            let res2 = fin2.fgets(&mut buffer2, MAX_LINE_SIZE);

            if !res && res2 {
                bomb_error("Paired input files have different number sequences");
            }
            if res && !res2 {
                bomb_error("Paired input files have different number sequences");
            }
            if one.size == 0 && two.size > 0 {
                bomb_error("Paired input files have different number sequences");
            }
            if one.size > 0 && two.size == 0 {
                bomb_error("Paired input files have different number sequences");
            }
            if !res && one.size == 0 {
                break;
            }

            if buffer[0] == b'+' {
                // fastq 3rd line
                // file 1
                let mut len = buffer.len();
                let mut len2 = len;
                while len2 != 0 && buffer[len2 - 1] != b'\n' {
                    // read until the end of the line
                    if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                        break;
                    }
                    len2 = buffer.len();
                    len += len2;
                }
                one.tot_length += len as i32;

                // read next line quality score
                if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                    bomb_error("can not read quality score after");
                }
                len = buffer.len();
                len2 = len;
                while len2 != 0 && buffer[len2 - 1] != b'\n' {
                    if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                        break;
                    }
                    len2 = buffer.len();
                    len += len2;
                }
                one.tot_length += len as i32;

                // file 2
                len = buffer2.len();
                len2 = len;
                while len2 != 0 && buffer2[len2 - 1] != b'\n' {
                    // read until the end of the line
                    if !fin2.fgets(&mut buffer2, MAX_LINE_SIZE) {
                        break;
                    }
                    len2 = buffer2.len();
                    len += len2;
                }
                two.tot_length += len as i32;

                // read next line quality score
                if !fin2.fgets(&mut buffer2, MAX_LINE_SIZE) {
                    bomb_error("can not read quality score after");
                }
                len = buffer2.len();
                len2 = len;
                while len2 != 0 && buffer2[len2 - 1] != b'\n' {
                    if !fin2.fgets(&mut buffer2, MAX_LINE_SIZE) {
                        break;
                    }
                    len2 = buffer2.len();
                    len += len2;
                }
                two.tot_length += len as i32;
            } else if buffer[0] == b'>' || buffer[0] == b'@' || (!res && one.size != 0) {
                if one.size != 0 && two.size != 0 {
                    // write previous record
                    if one.identifier.is_empty() || one.format() != 0 {
                        println!("Warning: from file \"{}\",", file);
                        println!(
                            "Discarding invalid sequence or sequence without identifier and description!\n"
                        );
                        if !one.identifier.is_empty() {
                            println!("{}", String::from_utf8_lossy(&one.identifier));
                        }
                        println!("{}", String::from_utf8_lossy(&one.data));
                        one.size = 0;
                        one.data.clear();
                        two.size = 0;
                        two.data.clear();
                    }
                    if two.identifier.is_empty() || two.format() != 0 {
                        println!("Warning: from file \"{}\",", file2);
                        println!(
                            "Discarding invalid sequence or sequence without identifier and description!\n"
                        );
                        if !two.identifier.is_empty() {
                            println!("{}", String::from_utf8_lossy(&two.identifier));
                        }
                        println!("{}", String::from_utf8_lossy(&two.data));
                        one.size = 0;
                        one.data.clear();
                        two.size = 0;
                        two.data.clear();
                    }
                    one.index = self.sequences.len() as i32;
                    if (one.size + two.size) > option_l {
                        if options.trim_len > 0 {
                            one.trim(options.trim_len);
                        }
                        if options.trim_len_r2 > 0 {
                            two.trim(options.trim_len_r2);
                        }
                        self.sequences.push(Sequence::from_pair(&one, &two, 1));
                    }
                }
                // R1
                one.size = 0;
                one.data.clear();
                one.tot_length = 0;

                let mut len = buffer.len();
                let mut len2 = len;
                des.size = 0;
                des.data.clear();
                des.append(&buffer);
                while len2 != 0 && buffer[len2 - 1] != b'\n' {
                    if !fin.fgets(&mut buffer, MAX_LINE_SIZE) {
                        break;
                    }
                    des.append(&buffer);
                    len2 = buffer.len();
                    len += len2;
                }
                let mut offset = fin.tell();
                one.des_begin = offset - len as u64; // offset of ">" or "@"
                one.tot_length += len as i32; // count first line

                let mut i = 0usize;
                if des.data_at(i) == b'>' || des.data_at(i) == b'@' || des.data_at(i) == b'+' {
                    i += 1;
                }
                if des.data_at(i) == b' ' || des.data_at(i) == b'\t' {
                    i += 1;
                }
                if options.des_len != 0 && options.des_len < des.size {
                    des.size = options.des_len;
                }
                while i < des.size as usize && !c_isspace(des.data_at(i)) {
                    i += 1;
                }
                // find first non-space letter
                one.identifier.clear();
                one.identifier
                    .extend_from_slice(&des.data[..i.min(des.data.len())]);

                // R2
                two.size = 0;
                two.data.clear();
                two.tot_length = 0;

                len = buffer2.len();
                len2 = len;
                while len2 != 0 && buffer2[len2 - 1] != b'\n' {
                    if !fin2.fgets(&mut buffer2, MAX_LINE_SIZE) {
                        break;
                    }
                    len2 = buffer2.len();
                    len += len2;
                }
                offset = fin2.tell();
                two.des_begin = offset - len as u64;
                two.tot_length += len as i32; // count first line
                two.identifier.clear();
                two.identifier.extend_from_slice(&one.identifier);
            } else {
                one.tot_length += buffer.len() as i32;
                one.append(&buffer);
                two.tot_length += buffer2.len() as i32;
                two.append(&buffer2);
            }
        }
    }

    /// `SequenceDB::WriteClusters` / `WriteClustersgz`
    ///
    /// COMPAT: the single-file writer ignores `-sf` (`sort_outputf`); only the paired-end
    /// writer honours it. Preserved.
    pub fn write_clusters(&self, db: &str, newdb: &str, _options: &Options) {
        let n = self.rep_seqs.len();
        let mut sorting: Vec<u64> = (0..n)
            .map(|i| {
                ((self.sequences[self.rep_seqs[i] as usize].index as u64) << 32)
                    | self.rep_seqs[i] as u64
            })
            .collect();
        sorting.sort_unstable();

        let mut fin = RecordSource::open(db);
        let fout = match std::fs::File::create(newdb) {
            Ok(f) => f,
            Err(_) => bomb_error("file opening failed"),
        };
        let mut fout = std::io::BufWriter::with_capacity(1 << 16, fout);

        let records: Vec<(u64, i32)> = sorting
            .iter()
            .map(|&s| {
                let seq = &self.sequences[(s & 0xffffffff) as usize];
                (seq.des_begin, seq.tot_length)
            })
            .collect();
        copy_records(&mut fin, &mut fout, &records);
    }

    /// `SequenceDB::WriteClusters( db, db_pe, newdb, newdb_pe, options )`
    pub fn write_clusters_pe(
        &self,
        db: &str,
        db_pe: &str,
        newdb: &str,
        newdb_pe: &str,
        options: &Options,
    ) {
        let n = self.rep_seqs.len();
        let mut sorting: Vec<u64> = (0..n)
            .map(|i| {
                ((self.sequences[self.rep_seqs[i] as usize].index as u64) << 32)
                    | self.rep_seqs[i] as u64
            })
            .collect();
        sorting.sort_unstable();

        //sort fasta / fastq
        let mut clstr_idx1: Vec<i32> = Vec::new();
        if options.sort_outputf != 0 {
            let mut clstr_size = vec![0i32; n];
            clstr_idx1 = (0..n as i32).collect();
            for seq in self.sequences.iter() {
                let id = seq.cluster_id;
                if id < 0 || id >= n as i32 {
                    continue;
                }
                clstr_size[id as usize] += 1;
            }
            quick_sort_idxr(&mut clstr_size, &mut clstr_idx1, 0, n as i32 - 1);
        }

        let pick = |i: usize| -> &Sequence {
            if options.sort_outputf != 0 {
                &self.sequences[self.rep_seqs[clstr_idx1[i] as usize] as usize]
            } else {
                &self.sequences[(sorting[i] & 0xffffffff) as usize]
            }
        };

        let mut fin = RecordSource::open(db);
        let mut fin_pe = RecordSource::open(db_pe);
        let fout = match std::fs::File::create(newdb) {
            Ok(f) => f,
            Err(_) => bomb_error("file opening failed"),
        };
        let fout_pe = match std::fs::File::create(newdb_pe) {
            Ok(f) => f,
            Err(_) => bomb_error("file opening failed"),
        };
        let mut fout = std::io::BufWriter::with_capacity(1 << 16, fout);
        let mut fout_pe = std::io::BufWriter::with_capacity(1 << 16, fout_pe);

        let r1: Vec<(u64, i32)> = (0..n).map(|i| (pick(i).des_begin, pick(i).tot_length)).collect();
        let r2: Vec<(u64, i32)> = (0..n)
            .map(|i| (pick(i).des_begin2, pick(i).tot_length2))
            .collect();
        copy_records(&mut fin, &mut fout, &r1);
        copy_records(&mut fin_pe, &mut fout_pe, &r2);
    }

    /// `SequenceDB::DivideSave`
    pub fn divide_save(&self, db: &str, newdb: &str, n: i32, _options: &Options) {
        if n == 0 || self.sequences.is_empty() {
            return;
        }

        let mut max_seg = self.total_letter as usize / n as usize + self.sequences[0].size as usize;
        if max_seg >= MAX_BIN_SWAP as usize {
            max_seg = MAX_BIN_SWAP as usize;
        }

        let mut fin = RecordSource::open(db);
        let mut seg_size: usize = 0;
        let mut seg = 0;
        let mut fout = std::io::BufWriter::with_capacity(
            1 << 16,
            std::fs::File::create(format!("{}-{}", newdb, 0))
                .unwrap_or_else(|_| bomb_error("file opening failed")),
        );
        let n = self.sequences.len();
        for i in 0..n {
            let seq = &self.sequences[i];

            seg_size += seq.size as usize;
            if seg_size >= max_seg {
                seg += 1;
                let _ = std::io::Write::flush(&mut fout);
                fout = std::io::BufWriter::with_capacity(
                    1 << 16,
                    std::fs::File::create(format!("{}-{}", newdb, seg))
                        .unwrap_or_else(|_| bomb_error("file opening failed")),
                );
                seg_size = seq.size as usize;
            }

            copy_records(&mut fin, &mut fout, &[(seq.des_begin, seq.tot_length)]);
        }
        let _ = std::io::Write::flush(&mut fout);
    }

    /// `SequenceDB::WriteExtra1D`
    pub fn write_extra_1d(&self, options: &Options) {
        let db_clstr = format!("{}.clstr", options.output);
        let db_clstr_bak = format!("{}.bak.clstr", options.output);
        let n = self.sequences.len();
        let mut sorting: Vec<i64> = (0..n)
            .map(|i| ((self.sequences[i].index as i64) << 32) | i as i64)
            .collect();
        sorting.sort_unstable();

        if options.backup_file {
            let f = std::fs::File::create(&db_clstr_bak).unwrap_or_else(|_| {
                bomb_error("file opening failed");
            });
            let mut fout = std::io::BufWriter::new(f);
            for &s in sorting.iter() {
                let seq = &self.sequences[(s & 0xffffffff) as usize];
                seq.print_info(seq.cluster_id, &mut fout, options);
            }
        }

        println!("writing clustering information");
        let m = self.rep_seqs.len();
        let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); m];
        for &s in sorting.iter() {
            let k = (s & 0xffffffff) as usize;
            let id = self.sequences[k].cluster_id;
            clusters[id as usize].push(k);
        }

        let f = std::fs::File::create(&db_clstr).unwrap_or_else(|_| {
            bomb_error("file opening failed");
        });
        let mut fout = std::io::BufWriter::with_capacity(1 << 16, f);

        if options.sort_output != 0 {
            let mut clstr_size: Vec<i32> = clusters.iter().map(|c| c.len() as i32).collect();
            let mut clstr_idx1: Vec<i32> = (0..m as i32).collect();
            quick_sort_idxr(&mut clstr_size, &mut clstr_idx1, 0, m as i32 - 1);

            for i in 0..m {
                let i0 = clstr_idx1[i] as usize;
                let _ = writeln!(fout, ">Cluster {}", i);
                for (k, &member) in clusters[i0].iter().enumerate() {
                    self.sequences[member].print_info(k as i32, &mut fout, options);
                }
            }
        } else {
            for i in 0..m {
                let _ = writeln!(fout, ">Cluster {}", i);
                for (k, &member) in clusters[i].iter().enumerate() {
                    self.sequences[member].print_info(k as i32, &mut fout, options);
                }
            }
        }
    }

    /// `SequenceDB::WriteExtra2D`
    ///
    /// `self` is db2 (the sequences being compared in), `other` is db1 (the reference whose
    /// members head the clusters).
    pub fn write_extra_2d(&self, other: &SequenceDB, options: &Options) {
        let db_clstr = format!("{}.clstr", options.output);
        let db_clstr_bak = format!("{}.bak.clstr", options.output);
        let n = other.sequences.len();
        let n2 = self.sequences.len();
        let mut sorting: Vec<i64> = (0..n)
            .map(|i| ((other.sequences[i].index as i64) << 32) | i as i64)
            .collect();
        sorting.sort_unstable();

        if options.backup_file {
            let f = std::fs::File::create(&db_clstr_bak).unwrap_or_else(|_| {
                bomb_error("file opening failed");
            });
            let mut fout = std::io::BufWriter::new(f);
            for &sv in sorting.iter() {
                let seq = &other.sequences[(sv & 0xffffffff) as usize];
                seq.print_info(seq.cluster_id, &mut fout, options);
            }
            for seq in self.sequences.iter() {
                if seq.state & IS_REDUNDANT != 0 {
                    seq.print_info(seq.cluster_id, &mut fout, options);
                }
            }
        }

        println!("writing clustering information");
        let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n2 {
            let id = self.sequences[i].cluster_id;
            if self.sequences[i].state & IS_REDUNDANT != 0 {
                clusters[id as usize].push(i);
            }
        }

        let f = std::fs::File::create(&db_clstr).unwrap_or_else(|_| {
            bomb_error("file opening failed");
        });
        let mut fout = std::io::BufWriter::with_capacity(1 << 16, f);
        for i in 0..n {
            let seq = &other.sequences[i];
            let _ = writeln!(fout, ">Cluster {}", i);
            seq.print_info(0, &mut fout, options);
            for (k, &member) in clusters[i].iter().enumerate() {
                self.sequences[member].print_info(k as i32 + 1, &mut fout, options);
            }
        }
    }

    /// `SequenceDB::SortDivide`
    pub fn sort_divide(&mut self, options: &mut Options, sort: bool, alphabet: &Alphabet) {
        let n = self.sequences.len();
        self.total_letter = 0;
        self.total_desc = 0;
        self.max_len = 0;
        self.min_len = usize::MAX;
        for seq in self.sequences.iter_mut() {
            let len = seq.size as usize;
            self.total_letter += len as i64;
            if len > self.max_len {
                self.max_len = len;
            }
            if len < self.min_len {
                self.min_len = len;
            }
            if seq.swap.is_none() {
                seq.convert_bases(alphabet);
            }
            self.total_desc += seq.identifier.len() as i64;
        }
        options.max_entries = self.max_len * MAX_TABLE_SEQ;
        // The `sizeof(INTs) <= 2` warning cannot trigger: INTs is 32-bit under LONG_SEQ.

        println!(
            "longest and shortest : {} and {}",
            self.max_len, self.min_len
        );
        println!("Total letters: {}", self.total_letter);
        // END change all the NR_seq to iseq

        self.len_n50 = (self.max_len + self.min_len) / 2; // will be properly set, if sort is true;
        if sort {
            // **************************** Form NR_idx[], Sort them from Long to short
            let mut sum: i64 = 0;
            let m = self.max_len - self.min_len + 1;
            let mut count = vec![0i32; m]; // count for each size = max_len - i
            let mut accum = vec![0i32; m]; // count for all size > max_len - i
            let mut offset = vec![0i32; m]; // offset from accum[i] when filling sorting
            let mut sorting: Vec<usize> = vec![0; n];

            for seq in self.sequences.iter() {
                count[self.max_len - seq.size as usize] += 1;
            }
            for i in 1..m {
                accum[i] = accum[i - 1] + count[i - 1];
            }
            for i in 0..m {
                sum += (self.max_len - i) as i64 * count[i] as i64;
                if sum >= (self.total_letter >> 1) {
                    self.len_n50 = self.max_len - i;
                    break;
                }
            }
            for i in 0..n {
                let len = self.max_len - self.sequences[i].size as usize;
                let id = accum[len] + offset[len];
                sorting[id as usize] = i;
                offset[len] += 1;
            }
            options.max_entries = 0;
            let mut permuted: Vec<Sequence> = Vec::with_capacity(n);
            for i in 0..n {
                permuted.push(std::mem::take(&mut self.sequences[sorting[i]]));
            }
            self.sequences = permuted;
            for i in 0..n {
                if i < MAX_TABLE_SEQ {
                    options.max_entries += self.sequences[i].size as usize;
                }
            }

            println!("Sequences have been sorted");
            // END sort them from long to short
        }
    }
}

//quick_sort_idx calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
/// `quick_sort_idx`
pub fn quick_sort_idx(a: &mut [i32], idx: &mut [i32], lo0: i32, hi0: i32) -> i32 {
    let mut lo = lo0;
    let mut hi = hi0;

    if hi0 > lo0 {
        let mid = a[((lo0 + hi0) / 2) as usize];

        while lo <= hi {
            while (lo < hi0) && (a[lo as usize] < mid) {
                lo += 1;
            }
            while (hi > lo0) && (a[hi as usize] > mid) {
                hi -= 1;
            }
            if lo <= hi {
                a.swap(lo as usize, hi as usize);
                idx.swap(lo as usize, hi as usize);
                lo += 1;
                hi -= 1;
            }
        } // while

        if lo0 < hi {
            quick_sort_idx(a, idx, lo0, hi);
        }
        if lo < hi0 {
            quick_sort_idx(a, idx, lo, hi0);
        }
    } // if ( hi0 > lo0)
    0
} // quick_sort_idx

//decreasing can not use reverse of quick_sort_idx due to tie
//quick_sort_idxr calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
/// `quick_sort_idxr`
pub fn quick_sort_idxr(a: &mut [i32], idx: &mut [i32], lo0: i32, hi0: i32) -> i32 {
    let mut lo = lo0;
    let mut hi = hi0;

    if hi0 > lo0 {
        let mid = a[((lo0 + hi0) / 2) as usize];

        while lo <= hi {
            while (lo < hi0) && (a[lo as usize] > mid) {
                lo += 1;
            }
            while (hi > lo0) && (a[hi as usize] < mid) {
                hi -= 1;
            }
            if lo <= hi {
                a.swap(lo as usize, hi as usize);
                idx.swap(lo as usize, hi as usize);
                lo += 1;
                hi -= 1;
            }
        } // while

        if lo0 < hi {
            quick_sort_idxr(a, idx, lo0, hi);
        }
        if lo < hi0 {
            quick_sort_idxr(a, idx, lo, hi0);
        }
    } // if ( hi0 > lo0)
    0
} // quick_sort_idxr

impl WorkingParam {
    /// `WorkingParam::WorkingParam` / `WorkingParam::Set`
    pub fn new(a1: f64, a2: f64, an: f64) -> WorkingParam {
        let mut p = WorkingParam::default();
        p.set(a1, a2, an);
        p
    }

    /// `WorkingParam::Set`
    pub fn set(&mut self, a1: f64, a2: f64, an: f64) {
        self.aa1_cutoff = a1;
        self.aas_cutoff = a2;
        self.aan_cutoff = an;
        self.len_upper_bound = 0;
        self.len_lower_bound = 0;
    }

    /// `WorkingParam::ControlShortCoverage`
    pub fn control_short_coverage(&mut self, len: i32, options: &Options) {
        self.len_eff = len;
        self.aln_cover_flag = 0;
        if (options.short_coverage > 0.0) || (options.min_control > 0) {
            // has alignment coverage control
            self.aln_cover_flag = 1;
            self.min_aln_len_s = (len as f64 * options.short_coverage) as i32;
            if len - options.short_control > self.min_aln_len_s {
                self.min_aln_len_s = len - options.short_control;
            }
            if options.min_control > self.min_aln_len_s {
                self.min_aln_len_s = options.min_control;
            }
        }
        if !options.global_identity {
            self.len_eff = self.min_aln_len_s; //global_identity==0
        }
    }

    /// `WorkingParam::ControlLongCoverage`
    pub fn control_long_coverage(&mut self, len2: i32, options: &Options) {
        if self.aln_cover_flag != 0 {
            self.min_aln_len_l = (len2 as f64 * options.long_coverage) as i32;
            if len2 - options.long_control > self.min_aln_len_l {
                self.min_aln_len_l = len2 - options.long_control;
            }
            if options.min_control > self.min_aln_len_l {
                self.min_aln_len_l = options.min_control;
            }
        }
    }

    /// `WorkingParam::ComputeRequiredBases`
    pub fn compute_required_bases(&mut self, naa: i32, ss: i32, options: &Options) {
        // d: distance, fraction of errors;
        // e: number of errors;
        // g: length of the maximum gap;
        // m: word length;
        // n: sequence length;
        // alignment length = n - g + 1;
        // d = e / (n - g + 1);
        // e >= 1, so that, g <= n + 1 - 1/d
        // word count = (n - g - m + 1) - (e - 1)*m;
        //            = (n - g - m + 1) - (d*(n - g + 1) - 1)*m
        //            = (n - g + 1) - d*m*(n - g + 1)
        //            = (n - g + 1)*(1 - d*m)
        // minimum word count is reached when g == n + 1 - 1/d
        // so, minimum word count = 1/d - m.
        // if g == band_width: word count = (n - band + 1)*(1 - d*m);
        if options.use_distance {
            let band = options.band_width + 1;
            let invd = (1.0 / (options.distance_thd + 1e-9)) as i32;
            let k = if self.len_eff < invd { self.len_eff } else { invd };
            let ks = self.len_eff - ss + 1;
            let kn = self.len_eff - naa + 1;
            let ks2 = invd - ss;
            let kn2 = invd - naa;
            let _ks3 =
                ((self.len_eff - band + 1) as f64 * (1.0 - options.distance_thd * ss as f64)) as i32;
            let _kn3 = ((self.len_eff - band + 1) as f64
                * (1.0 - options.distance_thd * naa as f64)) as i32;
            //if( ks3 > ks2 ) ks2 = ks3;
            //if( kn3 > kn2 ) kn2 = kn3;
            let _ = (ks, k);
            self.required_aa1 = if ks2 < k { ks2 } else { k };
            self.required_aas = self.required_aa1;
            self.required_aan = if kn2 < kn { kn2 } else { kn };
            if self.required_aa1 <= 0 {
                self.required_aa1 = 1;
                self.required_aas = 1;
            }
            if self.required_aan <= 0 {
                self.required_aan = 1;
            }
            //required_aa1 = required_aas = required_aan = 0;
            return;
        }
        // (N-K)-K*(1-C)*N = C*K*N-(K-1)*N-K = (C*K-K+1)*N-K
        self.required_aa1 = (self.len_eff - ss)
            - (ss as f64 * ((1.0 - self.aa1_cutoff) * self.len_eff as f64).ceil()) as i32;
        if self.required_aa1 < 0 {
            self.required_aa1 = 0;
        }
        self.required_aas = self.required_aa1;
        self.required_aan = (self.len_eff - naa)
            - (naa as f64 * ((1.0 - self.aa1_cutoff) * self.len_eff as f64).ceil()) as i32;
        if self.required_aan < 0 {
            self.required_aan = 0;
        }

        let aa1_old = (self.aa1_cutoff * self.len_eff as f64) as i32 - ss + 1;
        let aas_old = (self.aas_cutoff * self.len_eff as f64) as i32;
        let aan_old = (self.aan_cutoff * self.len_eff as f64) as i32;

        let thd = options.cluster_thd;
        //double rest = (len_eff - ss) / double(len_eff * ss);
        let rest = (self.len_eff - naa) as f64 / (self.len_eff * naa) as f64;
        let thd0 = 1.0 - rest;
        let mut fnew = 0.0;
        let mut fold = 1.0;
        if thd > thd0 {
            fnew = (thd - thd0) / rest;
            fold = 1.0 - fnew;
        }

        self.required_aa1 = (fnew * self.required_aa1 as f64 + fold * aa1_old as f64) as i32;
        self.required_aas = (fnew * self.required_aas as f64 + fold * aas_old as f64) as i32;
        self.required_aan = (fnew * self.required_aan as f64 + fold * aan_old as f64) as i32;
    }
}

// when alignment coverage such as -aL is specified
// if a existing rep is too long, it won't be qulified
/// `upper_bound_length_rep( int, double, int, double, int )`
pub fn upper_bound_length_rep_raw(len: i32, opt_s: f64, opt_s_aa: i32, opt_al: f64, opt_al_aa: i32) -> i32 {
    let mut len_upper_bound = 99999999;
    let r1 = if opt_s > opt_al { opt_s } else { opt_al };
    let a2 = if opt_s_aa < opt_al_aa { opt_s_aa } else { opt_al_aa };
    if r1 > 0.0 {
        len_upper_bound = (len as f32 / r1 as f32) as i32;
    }
    if (len + a2) < len_upper_bound {
        len_upper_bound = len + a2;
    }

    len_upper_bound
} // END upper_bound_length_rep

/// `upper_bound_length_rep( int, const Options & )`
pub fn upper_bound_length_rep(len: i32, options: &Options) -> i32 {
    upper_bound_length_rep_raw(
        len,
        options.diff_cutoff,
        options.diff_cutoff_aa,
        options.long_coverage,
        options.long_control,
    )
}

/// `cal_aax_cutoff`
pub fn cal_aax_cutoff(
    aa1_cutoff: &mut f64,
    aa2_cutoff: &mut f64,
    aan_cutoff: &mut f64,
    cluster_thd: f64,
    tolerance: i32,
    naa_stat_start_percent: i32,
    naa_stat: &[[[i32; 4]; 61]; 5],
    naa: i32,
) {
    *aa1_cutoff = cluster_thd;
    *aa2_cutoff = 1.0 - (1.0 - cluster_thd) * 2.0;
    *aan_cutoff = 1.0 - (1.0 - cluster_thd) * naa as f64;
    if tolerance == 0 {
        return;
    }

    let mut clstr_idx = (cluster_thd * 100.0) as i32 - naa_stat_start_percent;
    if clstr_idx < 0 {
        clstr_idx = 0;
    }
    let d2 = naa_stat[(tolerance - 1) as usize][clstr_idx as usize][3] as f64 / 100.0;
    let dn = naa_stat[(tolerance - 1) as usize][clstr_idx as usize][(5 - naa) as usize] as f64 / 100.0;
    *aa2_cutoff = if d2 > *aa2_cutoff { d2 } else { *aa2_cutoff };
    *aan_cutoff = if dn > *aan_cutoff { dn } else { *aan_cutoff };
} // END cal_aax_cutoff

/// `update_aax_cutoff`
pub fn update_aax_cutoff(
    aa1_cutoff: &mut f64,
    aa2_cutoff: &mut f64,
    aan_cutoff: &mut f64,
    tolerance: i32,
    naa_stat_start_percent: i32,
    naa_stat: &[[[i32; 4]; 61]; 5],
    naa: i32,
    cluster_thd: f64,
) {
    let cluster_thd = if cluster_thd > 1.0 { 1.00 } else { cluster_thd };

    let mut aa1_t = 0.0;
    let mut aa2_t = 0.0;
    let mut aan_t = 0.0;
    cal_aax_cutoff(
        &mut aa1_t,
        &mut aa2_t,
        &mut aan_t,
        cluster_thd,
        tolerance,
        naa_stat_start_percent,
        naa_stat,
        naa,
    );
    if aa1_t > *aa1_cutoff {
        *aa1_cutoff = aa1_t;
    }
    if aa2_t > *aa2_cutoff {
        *aa2_cutoff = aa2_t;
    }
    if aan_t > *aan_cutoff {
        *aan_cutoff = aan_t;
    }
} // END update_aax_cutoff

impl WorkingBuffer {
    /// `WorkingBuffer::WorkingBuffer`
    ///
    /// `maxlen` is the longest *query*, which is what upstream passes and what the
    /// `total_bytes` tally is based on. `max_diag_len` additionally bounds the longest
    /// *representative*: `diag_test_*` indexes its diagonal arrays by `len1 + len2 - 1`, so
    /// in 2D mode, where the representatives come from the other database, those arrays
    /// have to be sized from both. Upstream sidesteps this by fixing them at `MAX_DIAG`.
    pub fn new(
        frag: usize,
        maxlen: usize,
        max_diag_len: usize,
        options: &Options,
    ) -> WorkingBuffer {
        let mut b = WorkingBuffer {
            taap: Vec::new(),
            word_encodes: Vec::new(),
            word_encodes_backup: Vec::new(),
            word_encodes_no: Vec::new(),
            aap_list: Vec::new(),
            aap_begin: Vec::new(),
            look_counts: LookCounts::new(0),
            index_mapping: Vec::new(),
            score_mat: Vec::new(),
            back_mat: Vec::new(),
            dp_stride: 0,
            diag_score: Vec::new(),
            aan_list_comp: Vec::new(),
            seqi_comp: Vec::new(),
            total_bytes: 0,
        };
        b.set(frag, maxlen, max_diag_len, options);
        // `seqi_comp.resize( MAX_SEQ )` upstream; sized from the input here instead.
        b.seqi_comp = vec![0u8; maxlen + 1];
        b
    }

    /// `WorkingBuffer::Set`
    pub fn set(&mut self, frag: usize, maxlen: usize, max_diag_len: usize, options: &Options) {
        let est = options.is_est;
        let mut m = (MAX_UAA * MAX_UAA) as usize;
        let max_len = maxlen;
        let mut band = max_len * max_len;
        if est {
            m *= m;
        }
        if band > options.band_width as usize {
            band = options.band_width as usize;
        }
        self.taap.resize(m, 0);
        self.aap_list.resize(max_len, 0);
        self.aap_begin.resize(m, 0);
        self.word_encodes.resize(max_len, 0);
        self.word_encodes_no.resize(max_len, 0);
        self.word_encodes_backup.resize(max_len, 0);
        /* each table can not contain more than MAX_TABLE_SEQ representatives or fragments! */
        let frag = if frag > MAX_TABLE_SEQ { MAX_TABLE_SEQ } else { frag };
        self.look_counts = LookCounts::new(frag + 2);
        self.index_mapping.resize(frag + 2, 0);
        // `diag_score.resize( MAX_DIAG )` upstream, where MAX_DIAG == 2 * MAX_SEQ; the
        // requirement is only `len1 + len2 - 1`, so it is sized from the input here. The
        // *accounting* below still uses MAX_DIAG, because `total_bytes` feeds the table
        // size limits that decide where clustering blocks are cut.
        let diag = 2 * max_diag_len.max(max_len) + 2;
        self.diag_score.resize(diag, 0);
        self.aan_list_comp.resize(max_len, 0);

        const MAX_DIAG: usize = (MAX_SEQ as usize) << 1;
        let mut tb = max_len as i32;
        let add = |tb: &mut i32, v: usize| *tb = (*tb as usize).wrapping_add(v) as i32;
        add(&mut tb, m * 4); // taap
        add(&mut tb, max_len * 4); // word_encodes
        add(&mut tb, max_len * 4); // word_encodes_backup
        add(&mut tb, MAX_DIAG * 4); // diag_score
        add(&mut tb, MAX_DIAG * 4); // diag_score2
        add(&mut tb, max_len * 4); // aan_list_comp
        add(&mut tb, max_len * 4); // word_encodes_no, sizeof(INTs) == 4
        add(&mut tb, max_len * 4); // aap_list
        add(&mut tb, m * 4); // aap_begin
        add(&mut tb, (frag + 2) * 4); // indexMapping
        add(&mut tb, (frag + 2) * 8); // lookCounts, sizeof(IndexCount) == 8
        add(&mut tb, max_len * (band * 4 + 16)); // score rows, sizeof(VectorInt) == 16
        add(&mut tb, max_len * (band * 4 + 16)); // back rows, sizeof(VectorInt64) == 16
        self.total_bytes = tb;
    }

    /// `WorkingBuffer::EncodeWords`
    pub fn encode_words(&mut self, seq: &Sequence, naa: i32, est: bool, naa_tab: &Naa) -> i32 {
        let seqi = &seq.data;
        let len = seq.size;
        // check_word_encodes
        let aan_no = len - naa + 1;
        let mut skip = 0;
        for j in 0..aan_no as usize {
            let word = &seqi[j..];
            let mut encode = 0i32;
            let mut k1 = naa - 1;
            for k in 0..naa as usize {
                encode = encode
                    .wrapping_add((word[k] as i32).wrapping_mul(naa_tab.naan_array[k1 as usize]));
                k1 -= 1;
            }
            self.word_encodes[j] = encode;
            self.word_encodes_backup[j] = encode;
        }

        if est {
            for j in 0..len as usize {
                if seqi[j] >= 4 {
                    // here N is 4
                    let i0 = if j as i32 - naa + 1 > 0 {
                        j as i32 - naa + 1
                    } else {
                        0
                    };
                    let i1 = if (j as i32) < aan_no { j as i32 } else { aan_no - 1 };
                    for i in i0..=i1 {
                        self.word_encodes[i as usize] = -1;
                    }
                }
            }
            for j in 0..aan_no as usize {
                skip += (self.word_encodes[j] == -1) as i32;
            }
        }

        self.word_encodes[..aan_no as usize].sort_unstable();
        for j in 0..aan_no as usize {
            self.word_encodes_no[j] = 1;
        }
        let mut j = aan_no - 1;
        while j > 0 {
            if self.word_encodes[j as usize] == self.word_encodes[(j - 1) as usize] {
                self.word_encodes_no[(j - 1) as usize] += self.word_encodes_no[j as usize];
                self.word_encodes_no[j as usize] = 0;
            }
            j -= 1;
        }
        skip
        // END check_word_encodes
    }

    /// `WorkingBuffer::ComputeAAP`
    pub fn compute_aap(&mut self, seqi: &[u8], size: i32, naa_tab: &Naa) {
        let naa1 = naa_tab.naa1();
        let naa2 = naa_tab.naa2() as usize;
        let len1 = size - 1;
        for sk in 0..naa2 {
            self.taap[sk] = 0;
        }
        for j1 in 0..len1 as usize {
            let c22 = (seqi[j1] as i32) * naa1 + seqi[j1 + 1] as i32;
            self.taap[c22 as usize] += 1;
        }
        let mut mm = 0;
        for sk in 0..naa2 {
            self.aap_begin[sk] = mm;
            mm += self.taap[sk] as INTs;
            self.taap[sk] = 0;
        }
        for j1 in 0..len1 as usize {
            let c22 = ((seqi[j1] as i32) * naa1 + seqi[j1 + 1] as i32) as usize;
            let at = self.aap_begin[c22] as usize + self.taap[c22] as usize;
            self.aap_list[at] = j1 as INTs;
            self.taap[c22] += 1;
        }
    }

    /// `WorkingBuffer::ComputeAAP2`
    pub fn compute_aap2(&mut self, seqi: &[u8], size: i32, naa_tab: &Naa) {
        let naa1 = naa_tab.naa1();
        let naa2 = naa_tab.naa2();
        let naa3 = naa_tab.naa3();
        let naa4 = naa_tab.naa4() as usize;
        let len1 = size - 3;
        for sk in 0..naa4 {
            self.taap[sk] = 0;
        }
        for j1 in 0..len1 as usize {
            if seqi[j1] >= 4 || seqi[j1 + 1] >= 4 || seqi[j1 + 2] >= 4 || seqi[j1 + 3] >= 4 {
                continue; //skip N
            }
            let c22 = (seqi[j1] as i32) * naa3
                + (seqi[j1 + 1] as i32) * naa2
                + (seqi[j1 + 2] as i32) * naa1
                + seqi[j1 + 3] as i32;
            self.taap[c22 as usize] += 1;
        }
        let mut mm = 0;
        for sk in 0..naa4 {
            self.aap_begin[sk] = mm;
            mm += self.taap[sk] as INTs;
            self.taap[sk] = 0;
        }
        for j1 in 0..len1 as usize {
            if seqi[j1] >= 4 || seqi[j1 + 1] >= 4 || seqi[j1 + 2] >= 4 || seqi[j1 + 3] >= 4 {
                continue; //skip N
            }
            let c22 = ((seqi[j1] as i32) * naa3
                + (seqi[j1 + 1] as i32) * naa2
                + (seqi[j1 + 2] as i32) * naa1
                + seqi[j1 + 3] as i32) as usize;
            let at = self.aap_begin[c22] as usize + self.taap[c22] as usize;
            self.aap_list[at] = j1 as INTs;
            self.taap[c22] += 1;
        }
    }
}

/// The file-scope objects that upstream reads directly from inside the clustering
/// functions: `options`, `mat`, the `NAA*` word-size globals and `Comp_AAN_idx`.
#[derive(Clone, Copy)]
pub struct Globals<'a> {
    pub options: &'a Options,
    pub mat: &'a ScoreMatrix,
    pub naa: &'a Naa,
    /// `Comp_AAN_idx`: word code -> reverse-complement word code. Empty when `-r 0`.
    pub comp_aan_idx: &'a [i32],
}

/// Capacity of a `Vector<T>` after `n` calls to `Vector::Append`, which reserves
/// `n + n/5 + 1` whenever the next push would reach capacity.
///
/// `MinimalMemory` reports `sequences.capacity()`, and that figure feeds `mem_need`, which
/// decides where clustering blocks are cut -- so the growth schedule has to be reproduced
/// rather than approximated by `Vec`'s own doubling.
pub fn cpp_vector_capacity(n_appends: usize) -> usize {
    let mut cap = 0usize;
    for n in 0..n_appends {
        if n + 1 >= cap {
            cap = n + n / 5 + 1;
        }
    }
    cap
}

impl SequenceDB {
    /// `SequenceDB::ClusterOne`
    pub fn cluster_one(
        &mut self,
        id: i32,
        table: &mut WordTable,
        param: &mut WorkingParam,
        buffer: &mut WorkingBuffer,
        g: Globals,
    ) {
        let options = g.options;
        if self.sequences[id as usize].state & IS_REDUNDANT != 0 {
            return;
        }
        let frag_size = options.frag_size;
        let naa = options.naa;
        let len = self.sequences[id as usize].size;
        let len_bound = upper_bound_length_rep(len, options);
        param.len_upper_bound = len_bound;

        let mut seq = std::mem::take(&mut self.sequences[id as usize]);
        let flag = check_one(&mut seq, &self.sequences, table, param, buffer, g);

        if flag == 0 {
            if (seq.identity > 0.0) && options.cluster_best {
                // because of the -g option, this seq is similar to seqs in old SEGs
                seq.state |= IS_REDUNDANT;
                seq.clear();
            } else {
                // else add to NR90 db
                let aan_no = len - naa + 1;
                let size = self.rep_seqs.len() as i32;
                self.rep_seqs.push(id);
                seq.cluster_id = size;
                seq.identity = 0.0;
                seq.state |= IS_REP;
                if frag_size != 0 {
                    /* not used for EST */
                    let frg1 = (len - naa) / frag_size + 1;
                    let mut we = std::mem::take(&mut buffer.word_encodes_backup);
                    let mut wen = std::mem::take(&mut buffer.word_encodes_no);
                    table.add_word_counts_frag(aan_no, &mut we, &mut wen, frg1, frag_size);
                    buffer.word_encodes_backup = we;
                    buffer.word_encodes_no = wen;
                } else {
                    table.add_word_counts(
                        aan_no,
                        &buffer.word_encodes,
                        &buffer.word_encodes_no,
                        table.sequences.len() as i32,
                        options.is_est,
                    );
                }
                table.sequences.push(id as usize);
                if frag_size != 0 {
                    while (table.sequences.len() as i32) < table.frag_count {
                        table.sequences.push(id as usize);
                    }
                }
            }
        }
        self.sequences[id as usize] = seq;
        if (id + 1) % 1000 == 0 {
            let size = self.rep_seqs.len();
            print!(".");
            let _ = std::io::stdout().flush();
            if (id + 1) % 10000 == 0 {
                println!(
                    "\r..........{:9}  finished  {:9}  clusters",
                    id + 1,
                    size
                );
            }
        }
    }

    /// `SequenceDB::MinimalMemory`
    pub fn minimal_memory(
        &self,
        frag_no: i32,
        bsize: i32,
        t: i32,
        options: &Options,
        extra: usize,
        comp_aan_idx_len: usize,
    ) -> usize {
        let n = self.sequences.len();
        let f = if (frag_no as usize) < MAX_TABLE_SEQ {
            frag_no as usize
        } else {
            MAX_TABLE_SEQ
        };
        let mut mem_need: usize = 0;
        let mega: usize = 1000000;
        let table = if t > 1 { 2usize } else { 1usize };

        /// `sizeof(Sequence)` on x86-64, which `MinimalMemory` reports directly.
        const SIZEOF_SEQUENCE: usize = 112;
        /// `sizeof(NVector<IndexCount>)`
        const SIZEOF_NVECTOR: usize = 16;

        println!("\nApproximated minimal memory consumption:");
        let mut mem = n * SIZEOF_SEQUENCE + self.total_desc as usize + n + extra;
        if !options.store_disk {
            mem += self.total_letter as usize + n;
        }
        println!("{:<16}: {}M", "Sequence", mem / mega);
        mem_need += mem;

        mem = bsize as usize;
        println!(
            "{:<16}: {} X {}M = {}M",
            "Buffer",
            t,
            mem / mega,
            t as usize * mem / mega
        );
        mem_need += t as usize * mem;

        mem = f * (8 + 8) + self.naan as usize * SIZEOF_NVECTOR;
        println!(
            "{:<16}: {} X {}M = {}M",
            "Table",
            table,
            mem / mega,
            table * mem / mega
        );
        mem_need += table * mem;

        mem = cpp_vector_capacity(n) * 8 + n * 4;
        mem += comp_aan_idx_len * 4;
        println!("{:<16}: {}M", "Miscellaneous", mem / mega);
        mem_need += mem;

        println!("{:<16}: {}M\n", "Total", mem_need / mega);

        if options.max_memory != 0 && options.max_memory < mem_need + 50 * table {
            bomb_error(&format!(
                "not enough memory, please set -M option greater than {}\n",
                50 * table + mem_need / mega
            ));
        }
        mem_need
    }
}

/// `MemoryLimit`
///
/// COMPAT: the result is dead upstream -- both `DoClustering` overloads assign it to a
/// local that is never read. Retained for the 1:1 mapping.
pub fn memory_limit(mem_need: usize, options: &Options) -> usize {
    let mut mem_limit = (options.max_memory.wrapping_sub(mem_need)) / 8;

    if options.max_memory == 0 {
        mem_limit = options.max_entries;
        if mem_limit > MAX_TABLE_SIZE {
            mem_limit = MAX_TABLE_SIZE;
        }
    }
    mem_limit
}

impl Options {
    /// `Options::ComputeTableLimits`
    pub fn compute_table_limits(
        &mut self,
        _min_len: i32,
        max_len: i32,
        typical_len: i32,
        mem_need: usize,
    ) {
        //liwz Fri Jan 15 15:44:47 PST 2016
        //T=1 scale=1
        //T=2 scale=0.6035
        //T=4 scale=0.375
        //T=8 scale=0.2392
        //T=16 scale=0.1562
        //T=32 scale=0.104
        //T=64 scale=0.0703

        let scale = 0.5 / self.threads as f64 + 0.5 / (self.threads as f64).sqrt();
        self.max_sequences = (scale * MAX_TABLE_SEQ as f64) as usize;
        self.max_entries = (scale
            * (500i32.wrapping_mul(max_len))
                .wrapping_add(500000i32.wrapping_mul(typical_len))
                .wrapping_add(50000000) as f64) as usize;
        if self.max_memory != 0 {
            let frac = self.max_sequences as f64 / self.max_entries as f64;
            self.max_entries = (self.max_memory.wrapping_sub(mem_need)) / 8;
            self.max_sequences = (self.max_entries as f64 * frac) as usize;
            if self.max_sequences < MAX_TABLE_SEQ / 100 {
                self.max_sequences = MAX_TABLE_SEQ / 100;
            }
            if self.max_sequences > MAX_TABLE_SEQ {
                self.max_sequences = MAX_TABLE_SEQ;
            }
        }
        println!("Table limit with the given memory limit:");
        println!("Max number of representatives: {}", self.max_sequences);
        println!("Max number of word counting entries: {}\n", self.max_entries);
    }
}

/// `SequenceDB::CheckOne`
pub fn check_one(
    seq: &mut Sequence,
    db: &[Sequence],
    table: &WordTable,
    param: &mut WorkingParam,
    buf: &mut WorkingBuffer,
    g: Globals,
) -> i32 {
    let len = seq.size;
    param.len_upper_bound = upper_bound_length_rep(len, g.options);
    if g.options.is_est {
        return check_one_est(seq, db, table, param, buf, g);
    }
    check_one_aa(seq, db, table, param, buf, g)
}

/// `SequenceDB::CheckOneAA`
pub fn check_one_aa(
    seq: &mut Sequence,
    db: &[Sequence],
    table: &WordTable,
    param: &mut WorkingParam,
    buf: &mut WorkingBuffer,
    g: Globals,
) -> i32 {
    let options = g.options;
    // COMPAT: upstream copies the three cutoffs into locals and feeds them to
    // `update_aax_cutoff` at the bottom of the candidate loop, but then calls
    // `ComputeRequiredBases`, which reads `param`'s own fields. The locals are therefore
    // write-only and the `-g 1` cutoff tightening has no effect. Reproduced as written.
    let mut aa1_cutoff = param.aa1_cutoff;
    let mut aa2_cutoff = param.aas_cutoff;
    let mut aan_cutoff = param.aan_cutoff;

    let len = seq.size;
    let mut flag = 0;
    let frag_size = options.frag_size;

    let naa = options.naa;
    let s = table.sequences.len();
    let mut len_eff = len;

    if s != 0 {
        let mut min = db[table.sequences[s - 1]].size;
        if min < len {
            if len as f64 * options.diff_cutoff2 > min as f64 {
                min = (len as f64 * options.diff_cutoff2) as i32;
            }
            if (len - options.diff_cutoff_aa2) > min {
                min = len - options.diff_cutoff_aa2;
            }
            len_eff = min;
        }
    }

    //liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
    //longer seqeunce * -aL -band_width
    if s != 0 {
        let min = db[table.sequences[s - 1]].size;
        let min_red = (min as f64 * options.long_coverage - options.band_width as f64) as i32;
        if len < min_red {
            return 0; // return flag=0
        }
    }

    param.control_short_coverage(len_eff, options);
    param.compute_required_bases(options.naa, 2, options);

    buf.encode_words(seq, options.naa, false, g.naa);

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if options.min_control > len {
        return 0; // return flag=0
    }

    // lookup_aan
    let aan_no = len - options.naa + 1;
    {
        let look_counts = &mut buf.look_counts;
        let index_mapping = &mut buf.index_mapping;
        table.count_words(
            aan_no,
            &buf.word_encodes,
            &buf.word_encodes_no,
            look_counts,
            index_mapping,
            false,
            param.required_aan,
        );
    }

    // contained_in_old_lib()
    let len_upper_bound = param.len_upper_bound;
    let len_lower_bound = param.len_lower_bound;
    let mut band_left = 0;
    let mut band_right = 0;
    let mut best_score = 0;
    let mut best_sum = 0;
    let mut alnln = 0;
    let mut tiden_no = 0;
    let mut band_center = 0;
    let mut distance = 0.0f32;
    let mut talign_info = [0i32; 5];
    let frg2 = if frag_size != 0 {
        (len - naa + options.band_width) / frag_size + 1 + 1
    } else {
        0
    };
    let mut has_aa2 = 0;

    let mut ic_pos = 0usize;
    while buf.look_counts.items[ic_pos].count != 0 {
        let ic = buf.look_counts.items[ic_pos];
        if frag_size == 0 {
            buf.index_mapping[ic.index as usize] = 0;
            if ic.count < param.required_aan {
                ic_pos += 1;
                continue;
            }
        }

        let rep_idx = table.sequences[ic.index as usize];
        let len2 = db[rep_idx].size;
        if len2 > len_upper_bound {
            ic_pos += 1;
            continue;
        }
        if options.has_2d && len2 < len_lower_bound {
            ic_pos += 1;
            continue;
        }
        if frag_size != 0 {
            let ims = ic.index as usize;
            let k = (len2 - naa) / frag_size + 1;
            let mut sum = 0;
            for j1 in 0..frg2 {
                let im = buf.index_mapping[ims + j1 as usize];
                if im != 0 {
                    sum += buf.look_counts.items[(im - 1) as usize].count;
                }
            }
            let mut count = sum;
            for j1 in frg2..k {
                let im1 = buf.index_mapping[ims + j1 as usize];
                let im2 = buf.index_mapping[ims + (j1 - frg2) as usize];
                if im1 != 0 {
                    sum += buf.look_counts.items[(im1 - 1) as usize].count;
                }
                if im2 != 0 {
                    sum -= buf.look_counts.items[(im2 - 1) as usize].count;
                }
                if sum > count {
                    count = sum;
                }
            }
            if count < param.required_aan {
                ic_pos += 1;
                continue;
            }
        }

        param.control_long_coverage(len2, options);

        if has_aa2 == 0 {
            // calculate AAP array
            let seqi = std::mem::take(&mut seq.data);
            buf.compute_aap(&seqi, seq.size, g.naa);
            seq.data = seqi;
            has_aa2 = 1;
        }

        let band_width1 = if options.band_width < len + len2 - 2 {
            options.band_width
        } else {
            len + len2 - 2
        };
        diag_test_aapn(
            g.naa.naa1(),
            &db[rep_idx].data,
            len,
            len2,
            buf,
            &mut best_sum,
            band_width1,
            &mut band_left,
            &mut band_center,
            &mut band_right,
            param.required_aa1,
            options,
        );
        if best_sum < param.required_aas {
            ic_pos += 1;
            continue;
        }

        // The two arms of the upstream `if (options.print || aln_cover_flag)` are identical.
        let rc = local_band_align(
            &seq.data,
            &db[rep_idx].data,
            len,
            len2,
            g.mat,
            &mut best_score,
            &mut tiden_no,
            &mut alnln,
            &mut distance,
            Some(&mut talign_info),
            band_left,
            band_center,
            band_right,
            buf,
            options,
        );
        if rc == FAILED_FUNC {
            ic_pos += 1;
            continue;
        }
        if tiden_no < param.required_aa1 {
            ic_pos += 1;
            continue;
        }
        let mut lens = len;
        if options.has_2d && len > len2 {
            lens = len2;
        }
        let len_eff1 = if !options.global_identity {
            alnln
        } else {
            lens - talign_info[4]
        };
        let tiden_pc = tiden_no as f32 / len_eff1 as f32;
        if options.use_distance {
            if distance as f64 > options.distance_thd {
                ic_pos += 1;
                continue;
            }
            if distance >= seq.distance {
                ic_pos += 1;
                continue;
            } // existing distance
        } else {
            if (tiden_pc as f64) < options.cluster_thd {
                ic_pos += 1;
                continue;
            }
            if tiden_pc <= seq.identity {
                ic_pos += 1;
                continue;
            } // existing iden_no
        }
        if param.aln_cover_flag != 0 {
            if talign_info[3] - talign_info[2] + 1 < param.min_aln_len_l {
                ic_pos += 1;
                continue;
            }
            if talign_info[1] - talign_info[0] + 1 < param.min_aln_len_s {
                ic_pos += 1;
                continue;
            }
        }
        if options.has_2d {
            seq.state |= IS_REDUNDANT;
        }
        flag = 1;
        seq.identity = tiden_pc;
        seq.cluster_id = db[rep_idx].cluster_id;
        seq.distance = distance;
        seq.coverage[0] = talign_info[0] + 1;
        seq.coverage[1] = talign_info[1] + 1;
        seq.coverage[2] = talign_info[2] + 1;
        seq.coverage[3] = talign_info[3] + 1;
        if !options.cluster_best {
            break;
        }
        update_aax_cutoff(
            &mut aa1_cutoff,
            &mut aa2_cutoff,
            &mut aan_cutoff,
            options.tolerance,
            NAA_STAT_START_PERCENT,
            &NAA_STAT,
            naa,
            tiden_pc as f64,
        );
        param.compute_required_bases(options.naa, 2, options);
        ic_pos += 1;
    }
    if frag_size != 0 {
        ic_pos = 0;
    }
    while buf.look_counts.items[ic_pos].count != 0 {
        let idx = buf.look_counts.items[ic_pos].index as usize;
        buf.index_mapping[idx] = 0;
        ic_pos += 1;
    }
    buf.look_counts.size = 0;
    if flag == 1 {
        // if similar to old one delete it
        if !options.cluster_best {
            seq.clear();
            seq.state |= IS_REDUNDANT;
        }
    }
    flag
}

/// `SequenceDB::CheckOneEST`
pub fn check_one_est(
    seq: &mut Sequence,
    db: &[Sequence],
    table: &WordTable,
    param: &mut WorkingParam,
    buf: &mut WorkingBuffer,
    g: Globals,
) -> i32 {
    let options = g.options;
    let len = seq.size;
    let mut flag = 0;
    let s = table.sequences.len();
    let mut len_eff = len;
    if s != 0 {
        let mut min = db[table.sequences[s - 1]].size;
        if min < len {
            if len as f64 * options.diff_cutoff2 > min as f64 {
                min = (len as f64 * options.diff_cutoff2) as i32;
            }
            if (len - options.diff_cutoff_aa2) > min {
                min = len - options.diff_cutoff_aa2;
            }
            len_eff = min;
        }
    }

    //liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
    //longer seqeunce * -aL -band_width
    if s != 0 {
        let min = db[table.sequences[s - 1]].size;
        let min_red = (min as f64 * options.long_coverage - options.band_width as f64) as i32;
        if len < min_red {
            return 0; // return flag=0
        }
    }

    param.control_short_coverage(len_eff, options);
    param.compute_required_bases(options.naa, 4, options);
    let skip = buf.encode_words(seq, options.naa, true, g.naa);
    param.required_aan -= skip;
    param.required_aas -= skip;
    param.required_aa1 -= skip;
    if param.required_aan <= 0 {
        param.required_aan = 1;
    }
    if param.required_aas <= 0 {
        param.required_aas = 1;
    }
    if param.required_aa1 <= 0 {
        param.required_aa1 = 1;
    }

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if options.min_control > len {
        return 0; // return flag=0
    }

    let aan_no = len - options.naa + 1;

    // contained_in_old_lib()
    let len_upper_bound = param.len_upper_bound;
    let len_lower_bound = param.len_lower_bound;
    let mut band_left = 0;
    let mut band_right = 0;
    let mut best_score = 0;
    let mut best_sum = 0;
    let mut alnln = 0;
    let mut tiden_no = 0;
    let mut band_center = 0;
    let mut distance = 0.0f32;
    let mut talign_info = [0i32; 5];

    // `seqi` in the C++: the forward sequence for `comp == 0`, its reverse complement for
    // `comp == 1`. Both buffers are moved out of their owners so that `buf` and `seq` stay
    // independently mutable while a strand is borrowed.
    let seq_data: Vec<u8> = std::mem::take(&mut seq.data);
    let mut seqi_comp: Vec<u8> = std::mem::take(&mut buf.seqi_comp);

    for comp in 0..2 {
        if comp != 0 {
            for j0 in 0..aan_no as usize {
                let j = buf.word_encodes[j0];
                if j < 0 {
                    buf.aan_list_comp[j0] = j;
                } else {
                    buf.aan_list_comp[j0] = g.comp_aan_idx[j as usize];
                }
            }
            make_comp_iseq(len as usize, &mut seqi_comp, &seq_data);
        }
        let seqi: &[u8] = if comp == 0 { &seq_data } else { &seqi_comp };
        let mut has_aas = 0;

        {
            let look_counts = &mut buf.look_counts;
            let index_mapping = &mut buf.index_mapping;
            if comp != 0 {
                table.count_words(
                    aan_no,
                    &buf.aan_list_comp,
                    &buf.word_encodes_no,
                    look_counts,
                    index_mapping,
                    true,
                    param.required_aan,
                );
            } else {
                table.count_words(
                    aan_no,
                    &buf.word_encodes,
                    &buf.word_encodes_no,
                    look_counts,
                    index_mapping,
                    true,
                    param.required_aan,
                );
            }
        }

        let mut ic_pos = 0usize;
        while buf.look_counts.items[ic_pos].count != 0 {
            let ic = buf.look_counts.items[ic_pos];
            buf.index_mapping[ic.index as usize] = 0;
            if ic.count < param.required_aan {
                ic_pos += 1;
                continue;
            }
            let rep_idx = table.sequences[ic.index as usize];

            let len2 = db[rep_idx].size;
            if len2 > len_upper_bound {
                ic_pos += 1;
                continue;
            }
            if options.has_2d && len2 < len_lower_bound {
                ic_pos += 1;
                continue;
            }

            param.control_long_coverage(len2, options);

            if has_aas == 0 {
                // calculate AAP array
                buf.compute_aap2(seqi, seq.size, g.naa);
                has_aas = 1;
            }

            let band_width1 = if options.band_width < len + len2 - 2 {
                options.band_width
            } else {
                len + len2 - 2
            };
            diag_test_aapn_est(
                g.naa.naa1(),
                &db[rep_idx].data,
                len,
                len2,
                buf,
                &mut best_sum,
                band_width1,
                &mut band_left,
                &mut band_center,
                &mut band_right,
                param.required_aa1,
                options,
            );
            if best_sum < param.required_aas {
                ic_pos += 1;
                continue;
            }

            let rc = local_band_align(
                seqi,
                &db[rep_idx].data,
                len,
                len2,
                g.mat,
                &mut best_score,
                &mut tiden_no,
                &mut alnln,
                &mut distance,
                Some(&mut talign_info),
                band_left,
                band_center,
                band_right,
                buf,
                options,
            );
            if (options.print != 0 || param.aln_cover_flag != 0) && comp != 0 {
                //return overlap region
                talign_info[0] = len - talign_info[0] - 1;
                talign_info[1] = len - talign_info[1] - 1;
            }
            if rc == FAILED_FUNC {
                ic_pos += 1;
                continue;
            }
            if tiden_no < param.required_aa1 {
                ic_pos += 1;
                continue;
            }
            if options.is_454 {
                if talign_info[2] != talign_info[0] {
                    ic_pos += 1;
                    continue;
                } // same start
                if talign_info[0] > 1 {
                    ic_pos += 1;
                    continue;
                } // one mismatch allowed at beginning
                if (len - talign_info[1]) > 2 {
                    ic_pos += 1;
                    continue;
                } // one mismatch allowed at end
            }

            let mut lens = len;
            if options.has_2d && len > len2 {
                lens = len2;
            }
            let len_eff1 = if !options.global_identity {
                alnln
            } else {
                lens - talign_info[4]
            };
            let tiden_pc = tiden_no as f32 / len_eff1 as f32;
            if options.use_distance {
                if distance as f64 > options.distance_thd {
                    ic_pos += 1;
                    continue;
                }
                if options.cluster_best && distance >= seq.distance {
                    ic_pos += 1;
                    continue;
                } // existing distance
            } else {
                if (tiden_pc as f64) < options.cluster_thd {
                    ic_pos += 1;
                    continue;
                }
                if options.cluster_best && tiden_pc < seq.identity {
                    ic_pos += 1;
                    continue;
                } // existing iden_no
            }
            if param.aln_cover_flag != 0 {
                if talign_info[3] - talign_info[2] + 1 < param.min_aln_len_l {
                    ic_pos += 1;
                    continue;
                }
                if comp != 0 {
                    if talign_info[0] - talign_info[1] + 1 < param.min_aln_len_s {
                        ic_pos += 1;
                        continue;
                    }
                } else if talign_info[1] - talign_info[0] + 1 < param.min_aln_len_s {
                    ic_pos += 1;
                    continue;
                }
            }
            if options.cluster_best
                && (tiden_pc - seq.identity).abs() < 1e-9
                && db[rep_idx].cluster_id >= seq.cluster_id
            {
                ic_pos += 1;
                continue;
            }
            if !options.cluster_best && flag != 0 && db[rep_idx].cluster_id >= seq.cluster_id {
                ic_pos += 1;
                continue;
            }
            flag = if comp != 0 { -1 } else { 1 };
            seq.identity = tiden_pc;
            seq.distance = distance;
            seq.cluster_id = db[rep_idx].cluster_id;
            seq.coverage[0] = talign_info[0] + 1;
            seq.coverage[1] = talign_info[1] + 1;
            seq.coverage[2] = talign_info[2] + 1;
            seq.coverage[3] = talign_info[3] + 1;
            if !options.cluster_best {
                break;
            }
            ic_pos += 1;
        }
        while buf.look_counts.items[ic_pos].count != 0 {
            let idx = buf.look_counts.items[ic_pos].index as usize;
            buf.index_mapping[idx] = 0;
            ic_pos += 1;
        }
        buf.look_counts.size = 0;
        if options.option_r == 0 {
            break;
        }
    }
    buf.seqi_comp = seqi_comp;
    seq.data = seq_data;

    if flag == 1 || flag == -1 {
        // if similar to old one delete it
        if !options.cluster_best {
            seq.clear();
            seq.state |= IS_REDUNDANT;
        }
        if flag == -1 {
            seq.state |= IS_MINUS_STRAND;
        } else {
            seq.state &= !IS_MINUS_STRAND;
        }
    }
    flag
}

/// A `&mut [Sequence]` that hands out one element at a time to whichever worker claims
/// its index, used by the parallel screening phase of `do_clustering`.
struct SeqCells<'a> {
    ptr: *mut Sequence,
    len: usize,
    _marker: std::marker::PhantomData<&'a mut [Sequence]>,
}

// SAFETY: `Sequence` is `Send`, and `get` is only ever called with an index that the
// atomic cursor has handed to exactly one worker, so no two threads alias an element.
unsafe impl Send for SeqCells<'_> {}
unsafe impl Sync for SeqCells<'_> {}

impl<'a> SeqCells<'a> {
    fn new(slice: &'a mut [Sequence]) -> SeqCells<'a> {
        SeqCells {
            ptr: slice.as_mut_ptr(),
            len: slice.len(),
            _marker: std::marker::PhantomData,
        }
    }

    /// # Safety
    /// The caller must guarantee that `i < self.len` and that no other live borrow of
    /// element `i` exists. The atomic cursor in `screen_against_table` provides this: each
    /// index is returned by `fetch_add` exactly once.
    unsafe fn get(&self, i: usize) -> &'a mut Sequence {
        debug_assert!(i < self.len);
        &mut *self.ptr.add(i)
    }
}

/// `SequenceDB::ComputeDistance`
///
/// COMPAT: unreachable upstream -- its only call site sits inside a `#if 0` block at the
/// top of `DoClustering`. It is translated for completeness, including the bug on the
/// reverse-strand pass, where the copy loop is written `for(j=0;j<len;j++) comseq.data[i] =
/// seq->data[len-i-1];`: it indexes with the *outer* loop variable `i`, so it writes the
/// same slot `len` times and never actually reverses anything. Reproduced, except that the
/// write is skipped when `i` is past the end of the buffer -- upstream would corrupt memory
/// there, and there is no defined behaviour to match.
impl SequenceDB {
    pub fn compute_distance(
        &self,
        options: &Options,
        mat: &ScoreMatrix,
        naa_tab: &Naa,
        comp_aan_idx: &[i32],
    ) {
        let n = self.sequences.len();
        let mut best_score = 0;
        let mut best_sum = 0;
        let mut band_left = 0;
        let mut band_center = 0;
        let mut band_right = 0;
        let mut tiden_no = 0;
        let mut alnln = 0;
        let mut talign_info = [0i32; 5];
        let mut distance = 0.0f32;
        let mut buf = WorkingBuffer::new(n, self.max_len, self.max_len, options);
        let g = Globals {
            options,
            mat,
            naa: naa_tab,
            comp_aan_idx,
        };

        let mut dists: Vec<Vec<f32>> = vec![vec![0.0; n]; n];

        let mut comseq = self.sequences[0].clone();

        for i in 0..n {
            let seq = &self.sequences[i];
            let len = seq.size;
            buf.encode_words(seq, options.naa, false, naa_tab);
            buf.compute_aap2(&seq.data, seq.size, naa_tab);
            dists[i][i] = 0.0;
            if (i + 1) % 1000 == 0 {
                println!("{:9}", i + 1);
            }
            for j in 0..i {
                let rep = &self.sequences[j];
                let len2 = rep.size;
                let band_width1 = if options.band_width < len + len2 - 2 {
                    options.band_width
                } else {
                    len + len2 - 2
                };
                diag_test_aapn_est(
                    naa_tab.naa1(), &rep.data, len, len2, &mut buf, &mut best_sum,
                    band_width1, &mut band_left, &mut band_center, &mut band_right, 0, options,
                );
                local_band_align(
                    &seq.data, &rep.data, len, len2, mat, &mut best_score, &mut tiden_no,
                    &mut alnln, &mut distance, Some(&mut talign_info),
                    band_left, band_center, band_right, &mut buf, options,
                );
                dists[seq.index as usize][rep.index as usize] = distance;
                dists[rep.index as usize][seq.index as usize] = distance;
            }
            if options.option_r == 0 {
                break;
            }
            comseq.index = seq.index;
            comseq.size = len;
            for _j in 0..len {
                // COMPAT: the upstream index is `i`, not `j`. See the note above.
                if i < comseq.data.len() && (len - i as i32 - 1) >= 0 {
                    comseq.data[i] = seq.data[(len - i as i32 - 1) as usize];
                }
            }
            buf.encode_words(&comseq, options.naa, false, naa_tab);
            buf.compute_aap2(&comseq.data, seq.size, naa_tab);
            for j in 0..i {
                let rep = &self.sequences[j];
                let len2 = rep.size;
                let band_width1 = if options.band_width < len + len2 - 2 {
                    options.band_width
                } else {
                    len + len2 - 2
                };
                diag_test_aapn_est(
                    naa_tab.naa1(), &rep.data, len, len2, &mut buf, &mut best_sum,
                    band_width1, &mut band_left, &mut band_center, &mut band_right, 0, options,
                );
                local_band_align(
                    &comseq.data, &rep.data, len, len2, mat, &mut best_score, &mut tiden_no,
                    &mut alnln, &mut distance, Some(&mut talign_info),
                    band_left, band_center, band_right, &mut buf, options,
                );
                if distance < dists[seq.index as usize][rep.index as usize] {
                    dists[seq.index as usize][rep.index as usize] = distance;
                    dists[rep.index as usize][seq.index as usize] = distance;
                }
            }
        }
        let _ = g;
        let output = format!("{}.dist", options.output);
        let f = std::fs::File::create(&output)
            .unwrap_or_else(|_| bomb_error("file opening failed"));
        let mut fout = std::io::BufWriter::with_capacity(1 << 16, f);
        let _ = write!(fout, "1");
        for i in 1..n {
            let _ = write!(fout, "\t{}", i + 1);
        }
        let _ = writeln!(fout);
        for i in 0..n {
            let _ = write!(fout, "{}", fmt_g(dists[i][0]));
            for j in 1..n {
                let _ = write!(fout, "\t{}", fmt_g(dists[i][j]));
            }
            let _ = writeln!(fout);
        }
    }
}

/// `printf("%g", x)`: shortest of `%e` and `%f` with six significant digits, trailing zeros
/// removed.
fn fmt_g(x: f32) -> String {
    let x = x as f64;
    if x == 0.0 {
        return "0".to_string();
    }
    let exp = x.abs().log10().floor() as i32;
    if !(-4..6).contains(&exp) {
        let mut m = format!("{:.5e}", x);
        if let Some(pos) = m.find('e') {
            let (mant, e) = m.split_at(pos);
            let mant = mant.trim_end_matches('0').trim_end_matches('.');
            let ev: i32 = e[1..].parse().unwrap_or(0);
            m = format!("{}e{}{:02}", mant, if ev < 0 { "-" } else { "+" }, ev.abs());
        }
        m
    } else {
        let dec = (5 - exp).max(0) as usize;
        let m = format!("{:.*}", dec, x);
        if m.contains('.') {
            m.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            m
        }
    }
}

/// `SequenceDB::DoClustering`
impl SequenceDB {
    /// Screen `sequences[from..to]` against `table`, whose representatives all live below
    /// `from`, spreading the work over `threads` workers.
    ///
    /// Every sequence in the range is independent: `check_one` reads the table and the
    /// representatives, and writes only to the sequence it was handed. The result is
    /// therefore identical no matter how the range is distributed, which is what makes the
    /// parallel path reproduce the single-threaded output exactly.
    fn screen_against_table(
        &mut self,
        from: usize,
        to: usize,
        table: &WordTable,
        params: &mut [WorkingParam],
        buffers: &mut [WorkingBuffer],
        g: Globals,
        threads: usize,
    ) {
        if from >= to {
            return;
        }
        if threads <= 1 {
            for j in from..to {
                if self.sequences[j].state & IS_REDUNDANT != 0 {
                    continue;
                }
                let mut seq = std::mem::take(&mut self.sequences[j]);
                check_one(
                    &mut seq,
                    &self.sequences,
                    table,
                    &mut params[0],
                    &mut buffers[0],
                    g,
                );
                self.sequences[j] = seq;
            }
            return;
        }

        let (head, tail) = self.sequences.split_at_mut(from);
        let cells = SeqCells::new(&mut tail[..to - from]);
        let cursor = std::sync::atomic::AtomicUsize::new(0);
        let n = to - from;
        let head: &[Sequence] = head;

        std::thread::scope(|scope| {
            for (param, buf) in params.iter_mut().zip(buffers.iter_mut()).take(threads) {
                let cells = &cells;
                let cursor = &cursor;
                scope.spawn(move || loop {
                    let k = cursor.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if k >= n {
                        break;
                    }
                    // SAFETY: `k < n` and each `k` is produced by `fetch_add` exactly once,
                    // so this worker holds the only reference to element `k`.
                    let seq = unsafe { cells.get(k) };
                    if seq.state & IS_REDUNDANT != 0 {
                        continue;
                    }
                    check_one(seq, head, table, param, buf, g);
                });
            }
        });
    }

    /// `SequenceDB::DoClustering( const Options & options )` and
    /// `SequenceDB::DoClustering( int T, const Options & options )`.
    ///
    /// COMPAT: upstream has two drivers. The single-threaded one is translated faithfully
    /// below. The OpenMP one is *not* reproduced: it decides block boundaries from
    /// `may_stop` / `self_stop` / `stop`, plain shared `int`s written and read by several
    /// threads under nothing stronger than `#pragma omp flush`, so its output varies from
    /// run to run at the same `-T` and differs from `-T 1`. There is no fixed behaviour to
    /// match. Instead, `-T` here parallelises the phase that is provably order-independent
    /// -- screening the remaining sequences against a finished word table -- so the result
    /// is bit-identical to `-T 1` for every thread count.
    pub fn do_clustering(
        &mut self,
        options: &Options,
        mat: &ScoreMatrix,
        naa_tab: &Naa,
        comp_aan_idx: &[i32],
    ) {
        let naa = options.naa;
        let mut aa1_cutoff = options.cluster_thd;
        let mut aas_cutoff = 1.0 - (1.0 - options.cluster_thd) * 4.0;
        let mut aan_cutoff = 1.0 - (1.0 - options.cluster_thd) * options.naa as f64;
        let seq_no = self.sequences.len();
        let mut frag_no = seq_no as i32;
        let frag_size = options.frag_size;
        let t = options.threads.max(1) as usize;

        if frag_size != 0 {
            frag_no = 0;
            for i in 0..seq_no {
                frag_no += (self.sequences[i].size - naa) / frag_size + 1;
            }
        }

        if !options.is_est {
            cal_aax_cutoff(
                &mut aa1_cutoff,
                &mut aas_cutoff,
                &mut aan_cutoff,
                options.cluster_thd,
                options.tolerance,
                NAA_STAT_START_PERCENT,
                &NAA_STAT,
                naa,
            );
        }

        let mut params: Vec<WorkingParam> =
            (0..t).map(|_| WorkingParam::new(aa1_cutoff, aas_cutoff, aan_cutoff)).collect();
        let mut buffers: Vec<WorkingBuffer> = (0..t)
            .map(|_| WorkingBuffer::new(frag_no as usize, self.max_len, self.max_len, options))
            .collect();

        let mut word_table = WordTable::new(options.naa, self.naan);

        let mem_need = self.minimal_memory(
            frag_no,
            buffers[0].total_bytes,
            1,
            options,
            0,
            comp_aan_idx.len(),
        );
        let _mem_limit = memory_limit(mem_need, options);
        let mega: usize = 1000000;
        let n = self.sequences.len();

        let mut tabsize: usize = 0;

        // COMPAT: upstream sizes the tables from `threads`, so `-T` changes where blocks
        // are cut and hence the clustering. Pinning this to 1 is what lets `-T` be a pure
        // speed knob here.
        let mut opts = options.clone();
        opts.threads = 1;
        opts.compute_table_limits(
            self.min_len as i32,
            self.max_len as i32,
            self.len_n50 as i32,
            mem_need,
        );

        let g = Globals {
            options,
            mat,
            naa: naa_tab,
            comp_aan_idx,
        };

        let mut i = 0usize;
        while i < n {
            let redundancy = (self.rep_seqs.len() as f32 + 1.0) / (i as f32 + 1.0);
            let mut m = i;
            let mut sum: usize = 0;
            let max_items = opts.max_entries;
            let max_seqs = opts.max_sequences;
            let mut items: usize = 0;

            // find a block from i to m, so that this block can fit into a word table
            //     ...
            //  i  ++++++++++++++++++++++++++
            //     ++++++++++++++++++++
            //     ++++++++++++++++
            //  m  +++++++++++++
            //     ...
            while m < n && ((sum as f32 * redundancy) as usize) < max_seqs && items < max_items {
                let seq = &self.sequences[m];
                if seq.state & IS_REDUNDANT == 0 {
                    items += (seq.size as f32 * redundancy) as usize;
                    sum += 1;
                }
                m += 1;
            }
            println!("\rcomparing sequences from  {:9}  to  {:9}", i, m);
            let _ = std::io::stdout().flush();
            for ks in i..m {
                // clustering this block
                i = ks + 1;
                if self.sequences[ks].state & IS_REDUNDANT != 0 {
                    continue;
                }
                self.cluster_one(ks as i32, &mut word_table, &mut params[0], &mut buffers[0], g);
                if word_table.size >= max_items {
                    break;
                }
                let tmax = max_seqs
                    - if frag_size != 0 {
                        (self.sequences[ks].size / frag_size + 1) as usize
                    } else {
                        0
                    };
                if word_table.sequences.len() >= tmax {
                    break;
                }
            } // finishing word table from this block
            m = i;
            if word_table.size == 0 {
                continue;
            }

            // use this word table to screen rest sequences m->N
            //
            // The serial loop stops at the first non-redundant `j` whose length upper bound
            // has fallen below the shortest representative in the table. Sequences are
            // sorted longest-first, so that bound only shrinks as `j` grows and the stopping
            // point is a property of the input, not of the iteration -- compute it up front
            // so the rest of the range can be screened in any order.
            let shortest_rep =
                self.sequences[word_table.sequences[word_table.sequences.len() - 1]].size;
            let mut j_end = n;
            for j in m..n {
                if self.sequences[j].state & IS_REDUNDANT != 0 {
                    continue;
                }
                if shortest_rep > upper_bound_length_rep(self.sequences[j].size, options) {
                    j_end = j + 1;
                    break;
                }
            }
            self.screen_against_table(m, j_end, &word_table, &mut params, &mut buffers, g, t);

            if word_table.size > tabsize {
                tabsize = word_table.size;
            }
            word_table.clear();
        }
        println!(
            "\n{:9}  finished  {:9}  clusters",
            self.sequences.len(),
            self.rep_seqs.len()
        );
        let mem = (mem_need + tabsize * 8) / mega;
        println!("\nApproximated maximum memory consumption: {}M", mem);
        word_table.clear();
    }
}

/// `SequenceDB::ClusterTo`
///
/// The 2D mode: `self` is db2, `other` is db1. db1's sequences become the representatives,
/// and every db2 sequence that matches one is marked redundant. db1 is never clustered.
impl SequenceDB {
    pub fn cluster_to(
        &mut self,
        other: &mut SequenceDB,
        options: &Options,
        mat: &ScoreMatrix,
        naa_tab: &Naa,
        comp_aan_idx: &[i32],
    ) {
        let naa = options.naa;
        let mut aa1_cutoff = options.cluster_thd;
        let mut aas_cutoff = 1.0 - (1.0 - options.cluster_thd) * 4.0;
        let mut aan_cutoff = 1.0 - (1.0 - options.cluster_thd) * options.naa as f64;
        let mut word_encodes: Vec<i32> = vec![0; other.max_len + 1];
        let mut word_encodes_no: Vec<INTs> = vec![0; other.max_len + 1];

        if !options.is_est {
            cal_aax_cutoff(
                &mut aa1_cutoff,
                &mut aas_cutoff,
                &mut aan_cutoff,
                options.cluster_thd,
                options.tolerance,
                NAA_STAT_START_PERCENT,
                &NAA_STAT,
                naa,
            );
        }

        let n = other.sequences.len();
        let m_total = self.sequences.len();
        let t = options.threads.max(1) as usize;
        let mut nr2_red_no = 0i32;

        let mut params: Vec<WorkingParam> = (0..t)
            .map(|_| WorkingParam::new(aa1_cutoff, aas_cutoff, aan_cutoff))
            .collect();
        let mut buffers: Vec<WorkingBuffer> = (0..t)
            .map(|_| WorkingBuffer::new(n, self.max_len, self.max_len.max(other.max_len), options))
            .collect();

        let mem_need = self.minimal_memory(
            n as i32,
            buffers[0].total_bytes,
            1,
            options,
            other.total_letter as usize + other.total_desc as usize,
            comp_aan_idx.len(),
        );
        let _mem_limit = memory_limit(mem_need, options);

        // COMPAT: as in `do_clustering`, the table limits are computed with `threads = 1`
        // so that `-T` cannot move the block boundaries and change the clustering.
        let mut opts = options.clone();
        opts.threads = 1;
        opts.compute_table_limits(
            self.min_len as i32,
            self.max_len as i32,
            self.len_n50 as i32,
            mem_need,
        );

        let g = Globals {
            options,
            mat,
            naa: naa_tab,
            comp_aan_idx,
        };

        let mut word_table = WordTable::new(options.naa, self.naan);

        let max_items = opts.max_entries;
        let max_seqs = opts.max_sequences;
        let mut i = 0usize;
        while i < n {
            let mut items: usize = 0;
            let mut sum: usize = 0;
            let mut m = i;
            while m < n && sum < max_seqs && items < max_items {
                let seq = &other.sequences[m];
                if seq.state & IS_REDUNDANT == 0 {
                    items += seq.size as usize;
                    sum += 1;
                }
                m += 1;
            }
            for ks in i..m {
                let seq = &mut other.sequences[ks];
                let len = seq.size;
                let mut aan_no = 0;
                calc_ann_list(
                    len,
                    &seq.data,
                    naa,
                    &mut aan_no,
                    &mut word_encodes,
                    &mut word_encodes_no,
                    options.is_est,
                    naa_tab,
                );
                seq.cluster_id = ks as i32;
                seq.state |= IS_REP;
                word_table.add_word_counts(
                    aan_no,
                    &word_encodes,
                    &word_encodes_no,
                    (ks - i) as i32,
                    options.is_est,
                );
                word_table.sequences.push(ks);
                if (ks + 1) % 1000 == 0 {
                    print!(".");
                    let _ = std::io::stdout().flush();
                    if (ks + 1) % 10000 == 0 {
                        println!("{:9}  finished", ks + 1);
                    }
                }
            }

            // The serial loop stops at the first non-redundant `j` whose length upper bound
            // has dropped below the shortest representative in the table. db2 is *not*
            // length-sorted here (`SortDivide` is called with `sort = false`), so the
            // stopping point is not monotone -- but it depends only on lengths and on the
            // redundancy flags as they stand now, neither of which the screening changes.
            // Find it first, then the range before it can be screened in any order.
            let shortest_rep =
                other.sequences[word_table.sequences[word_table.sequences.len() - 1]].size;
            let mut j_end = m_total;
            for j in 0..m_total {
                if self.sequences[j].state & IS_REDUNDANT != 0 {
                    continue;
                }
                if shortest_rep > upper_bound_length_rep(self.sequences[j].size, options) {
                    j_end = j;
                    break;
                }
            }

            nr2_red_no += screen_2d(
                &mut self.sequences[..j_end],
                &other.sequences,
                &word_table,
                &mut params,
                &mut buffers,
                g,
                t,
            );

            println!(
                "\r..........{:9}  compared  {:9}  clusters",
                i, nr2_red_no
            );
            word_table.clear();
            word_table.size = 0;
            i = m;
        }

        if options.cluster_best {
            //delete redundant sequences in options.cluster_best mode
            for seq in self.sequences.iter_mut() {
                if seq.identity > 0.0 {
                    seq.state |= IS_REDUNDANT;
                    nr2_red_no += 1;
                }
            }
        }
        for i in 0..self.sequences.len() {
            let seq = &mut self.sequences[i];
            if seq.identity < 0.0 {
                seq.identity *= -1.0;
            }
            if seq.state & IS_REDUNDANT == 0 {
                self.rep_seqs.push(i as i32);
            }
        }

        println!();
        println!(
            "{} compared\t{} clustered",
            self.sequences.len(),
            nr2_red_no
        );
    }
}

/// The screening half of `ClusterTo`: compare every db2 sequence in `queries` against the
/// db1 representatives in `table`, and return how many were found redundant.
///
/// Split out so the loop body is written once for the serial and threaded paths. As in
/// `do_clustering`, every query is independent -- it reads the finished table and writes
/// only to itself -- so the outcome does not depend on the number of workers.
fn screen_2d(
    queries: &mut [Sequence],
    db: &[Sequence],
    table: &WordTable,
    params: &mut [WorkingParam],
    buffers: &mut [WorkingBuffer],
    g: Globals,
    threads: usize,
) -> i32 {
    let options = g.options;
    let n = queries.len();

    let one = |seq: &mut Sequence, param: &mut WorkingParam, buf: &mut WorkingBuffer| -> i32 {
        let len = seq.size;
        let len_upper_bound = upper_bound_length_rep(len, options);
        let mut len_lower_bound = len - options.diff_cutoff_aa2;

        let len_tmp = (len as f64 * options.diff_cutoff2) as i32;
        if len_tmp < len_lower_bound {
            len_lower_bound = len_tmp;
        }
        param.len_upper_bound = len_upper_bound;
        param.len_lower_bound = len_lower_bound;

        let flag = check_one(seq, db, table, param, buf, g);
        let mut red = 0;
        if flag == 1 || flag == -1 {
            // if similar to old one delete it
            if !options.cluster_best {
                seq.clear();
                seq.state |= IS_REDUNDANT;
                red = 1;
            }
            if flag == -1 {
                seq.state |= IS_MINUS_STRAND; // for EST only
            }
        }
        red
    };

    if threads <= 1 || n < 64 {
        let (param, rest) = params.split_first_mut().unwrap();
        let _ = rest;
        let buf = &mut buffers[0];
        let mut red = 0;
        let mut p0 = 0.0f32;
        for (j, seq) in queries.iter_mut().enumerate() {
            if seq.state & IS_REDUNDANT != 0 {
                continue;
            }
            red += one(seq, param, buf);
            let p = (100.0 * j as f32) / n as f32;
            if p > p0 + 1e-1 {
                // print only if the percentage changed
                print!("\r{:4.1}%", p);
                let _ = std::io::stdout().flush();
                p0 = p;
            }
        }
        return red;
    }

    let cells = SeqCells::new(queries);
    let cursor = std::sync::atomic::AtomicUsize::new(0);
    let total = std::sync::atomic::AtomicI32::new(0);
    std::thread::scope(|scope| {
        for (param, buf) in params.iter_mut().zip(buffers.iter_mut()).take(threads) {
            let cells = &cells;
            let cursor = &cursor;
            let total = &total;
            let one = &one;
            scope.spawn(move || {
                let mut red = 0;
                loop {
                    let k = cursor.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if k >= n {
                        break;
                    }
                    // SAFETY: `k < n`, and `fetch_add` yields each `k` to one worker only.
                    let seq = unsafe { cells.get(k) };
                    if seq.state & IS_REDUNDANT != 0 {
                        continue;
                    }
                    red += one(seq, param, buf);
                }
                total.fetch_add(red, std::sync::atomic::Ordering::Relaxed);
            });
        }
    });
    total.load(std::sync::atomic::Ordering::Relaxed)
}

/// `calc_ann_list`
pub fn calc_ann_list(
    len: i32,
    seqi: &[u8],
    naa: i32,
    aan_no: &mut i32,
    aan_list: &mut [i32],
    aan_list_no: &mut [INTs],
    est: bool,
    naa_tab: &Naa,
) -> i32 {
    // check_word_encodes
    *aan_no = len - naa + 1;
    for j in 0..*aan_no as usize {
        let word = &seqi[j..];
        let mut encode = 0i32;
        let mut k1 = naa - 1;
        for k in 0..naa as usize {
            encode = encode
                .wrapping_add((word[k] as i32).wrapping_mul(naa_tab.naan_array[k1 as usize]));
            k1 -= 1;
        }
        aan_list[j] = encode;
    }
    if est {
        for j in 0..len as usize {
            if seqi[j] >= 4 {
                // here N is 4
                let i0 = if j as i32 - naa + 1 > 0 {
                    j as i32 - naa + 1
                } else {
                    0
                };
                let i1 = if (j as i32) < *aan_no {
                    j as i32
                } else {
                    *aan_no - 1
                };
                for i in i0..=i1 {
                    aan_list[i as usize] = -1;
                }
            }
        }
    }

    aan_list[..*aan_no as usize].sort_unstable();
    for j in 0..*aan_no as usize {
        aan_list_no[j] = 1;
    }
    let mut j = *aan_no - 1;
    while j > 0 {
        if aan_list[j as usize] == aan_list[(j - 1) as usize] {
            aan_list_no[(j - 1) as usize] += aan_list_no[j as usize];
            aan_list_no[j as usize] = 0;
        }
        j -= 1;
    }
    OK_FUNC
}

/// `make_comp_short_word_index`
pub fn make_comp_short_word_index(naa: i32, naan_array: &[i32; 13], comp_aan_idx: &mut [i32]) {
    let c: [usize; 4] = [3, 2, 1, 0];
    let mut short_word = [0u8; 32]; //short_word[12] is enough

    let naa1 = naan_array[1];
    let naan = naan_array[naa as usize];

    for i in 0..naan {
        // decompose i back to short_word
        let mut k = i;
        for j in 0..naa as usize {
            short_word[j] = (k % naa1) as u8;
            k /= naa1;
        }

        // calc_comp_aan_list
        let mut icomp = 0i32;
        let mut k1 = naa - 1;
        for k in 0..naa as usize {
            icomp += c[short_word[k1 as usize] as usize] as i32 * naan_array[k];
            k1 -= 1;
        }

        comp_aan_idx[i as usize] = icomp;
    }
} // make_comp_short_word_index

/// `naa_stat_start_percent`
pub const NAA_STAT_START_PERCENT: i32 = 40;

/// `naa_stat[tolerance-1][identity_pct - 40][5 - word_len]`: the empirically derived
/// minimum shared-word percentage for a given tolerance, identity and word length.
/// Columns are word lengths N=5, N=4, N=3, N=2.
pub static NAA_STAT: [[[i32; 4]; 61]; 5] = [
    // cover 0.99
    [
        //  N=5   N=4   N=3   N=2
        [   0,    0,    0,    7], // 40%
        [   0,    0,    0,    8], // 41%
        [   0,    0,    0,    9], // 42%
        [   0,    0,    0,    9], // 43%
        [   0,    0,    1,   10], // 44%
        [   0,    0,    1,   11], // 45%
        [   0,    0,    1,   12], // 46%
        [   0,    0,    2,   13], // 47%
        [   0,    0,    2,   14], // 48%
        [   0,    0,    4,   16], // 49%
        [   0,    0,    4,   16], // 50%
        [   0,    0,    5,   17], // 51%
        [   0,    0,    5,   18], // 52%
        [   0,    0,    7,   20], // 53%
        [   0,    1,    7,   21], // 54%
        [   0,    1,    7,   21], // 55%
        [   0,    2,    8,   23], // 56%
        [   0,    2,    8,   25], // 57%
        [   0,    2,   10,   25], // 58%
        [   0,    3,   10,   26], // 59%
        [   0,    4,   13,   28], // 60%
        [   0,    5,   13,   30], // 61%
        [   0,    5,   14,   30], // 62%
        [   1,    6,   15,   33], // 63%
        [   2,    7,   17,   34], // 64%
        [   2,    7,   17,   35], // 65%
        [   2,    9,   20,   37], // 66%
        [   4,   10,   20,   37], // 67%
        [   4,   11,   22,   40], // 68%
        [   5,   12,   24,   41], // 69%
        [   5,   12,   25,   42], // 70%
        [   6,   16,   27,   43], // 71%
        [   8,   16,   27,   45], // 72%
        [   9,   17,   29,   47], // 73%
        [  10,   18,   31,   47], // 74%
        [  10,   20,   32,   50], // 75%
        [  12,   20,   32,   51], // 76%
        [  14,   22,   36,   54], // 77%
        [  15,   24,   37,   55], // 78%
        [  17,   26,   41,   58], // 79%
        [  18,   29,   41,   59], // 80%
        [  20,   30,   45,   60], // 81%
        [  24,   35,   48,   62], // 82%
        [  26,   36,   48,   64], // 83%
        [  27,   38,   51,   65], // 84%
        [  31,   43,   54,   68], // 85%
        [  35,   43,   55,   70], // 86%
        [  36,   48,   60,   71], // 87%
        [  36,   50,   61,   73], // 88%
        [  40,   50,   61,   75], // 89%
        [  45,   54,   65,   75], // 90%
        [  52,   60,   70,   79], // 91%
        [  53,   62,   71,   81], // 92%
        [  57,   66,   75,   84], // 93%
        [  57,   66,   76,   85], // 94%
        [  64,   71,   78,   85], // 95%
        [  70,   75,   82,   89], // 96%
        [  77,   81,   86,   92], // 97%
        [  82,   86,   90,   94], // 98%
        [  83,   87,   91,   95], // 99%
        [  91,   93,   95,   97], // 100%
    ],
    // cover 0.95
    [
        //  N=5   N=4   N=3   N=2
        [   0,    0,    1,    9], // 40%
        [   0,    0,    2,   10], // 41%
        [   0,    0,    2,   11], // 42%
        [   0,    0,    3,   12], // 43%
        [   0,    0,    3,   12], // 44%
        [   0,    0,    4,   14], // 45%
        [   0,    0,    4,   14], // 46%
        [   0,    1,    5,   16], // 47%
        [   0,    1,    6,   17], // 48%
        [   0,    2,    7,   19], // 49%
        [   0,    2,    8,   19], // 50%
        [   0,    2,    8,   20], // 51%
        [   0,    2,    9,   21], // 52%
        [   0,    4,   10,   23], // 53%
        [   1,    4,   11,   24], // 54%
        [   1,    4,   11,   24], // 55%
        [   1,    5,   13,   26], // 56%
        [   2,    5,   13,   27], // 57%
        [   2,    6,   15,   29], // 58%
        [   2,    7,   15,   30], // 59%
        [   3,    8,   16,   31], // 60%
        [   4,    8,   18,   32], // 61%
        [   4,    9,   18,   33], // 62%
        [   5,   11,   20,   36], // 63%
        [   6,   12,   22,   37], // 64%
        [   6,   12,   22,   38], // 65%
        [   8,   14,   24,   40], // 66%
        [   8,   15,   25,   41], // 67%
        [  10,   16,   27,   42], // 68%
        [  10,   18,   28,   45], // 69%
        [  11,   18,   29,   45], // 70%
        [  14,   21,   31,   47], // 71%
        [  14,   22,   32,   48], // 72%
        [  14,   22,   33,   50], // 73%
        [  17,   24,   36,   52], // 74%
        [  17,   25,   36,   52], // 75%
        [  18,   27,   39,   54], // 76%
        [  20,   29,   41,   56], // 77%
        [  21,   31,   42,   58], // 78%
        [  21,   31,   46,   60], // 79%
        [  27,   35,   46,   60], // 80%
        [  28,   37,   50,   63], // 81%
        [  31,   38,   50,   64], // 82%
        [  34,   43,   53,   66], // 83%
        [  36,   45,   54,   67], // 84%
        [  41,   50,   60,   70], // 85%
        [  43,   51,   60,   71], // 86%
        [  45,   54,   63,   74], // 87%
        [  48,   55,   64,   75], // 88%
        [  54,   60,   68,   78], // 89%
        [  55,   62,   71,   80], // 90%
        [  56,   63,   71,   80], // 91%
        [  64,   70,   76,   84], // 92%
        [  69,   74,   80,   86], // 93%
        [  73,   78,   83,   88], // 94%
        [  74,   78,   84,   89], // 95%
        [  80,   84,   87,   91], // 96%
        [  83,   86,   90,   93], // 97%
        [  86,   89,   92,   95], // 98%
        [  91,   93,   95,   97], // 99%
        [  92,   93,   95,   97], // 100%
    ],
    // cover 0.90
    [
        //  N=5   N=4   N=3   N=2
        [   0,    0,    2,   11], // 40%
        [   0,    0,    3,   12], // 41%
        [   0,    0,    3,   12], // 42%
        [   0,    1,    4,   13], // 43%
        [   0,    1,    5,   14], // 44%
        [   0,    1,    5,   15], // 45%
        [   0,    1,    6,   16], // 46%
        [   0,    2,    7,   18], // 47%
        [   0,    2,    7,   18], // 48%
        [   0,    3,    9,   20], // 49%
        [   1,    4,    9,   20], // 50%
        [   1,    4,   10,   21], // 51%
        [   1,    4,   11,   23], // 52%
        [   2,    5,   12,   24], // 53%
        [   2,    5,   12,   25], // 54%
        [   2,    6,   13,   26], // 55%
        [   3,    7,   14,   28], // 56%
        [   3,    7,   15,   28], // 57%
        [   4,    8,   16,   30], // 58%
        [   5,    9,   17,   31], // 59%
        [   5,   10,   18,   32], // 60%
        [   6,   11,   20,   35], // 61%
        [   6,   11,   20,   35], // 62%
        [   7,   13,   22,   38], // 63%
        [   8,   14,   23,   39], // 64%
        [   8,   15,   24,   39], // 65%
        [  10,   16,   26,   42], // 66%
        [  10,   17,   27,   42], // 67%
        [  12,   19,   29,   44], // 68%
        [  13,   20,   30,   46], // 69%
        [  13,   21,   31,   47], // 70%
        [  16,   23,   33,   48], // 71%
        [  18,   25,   34,   50], // 72%
        [  18,   26,   36,   51], // 73%
        [  19,   28,   38,   53], // 74%
        [  20,   29,   38,   53], // 75%
        [  23,   30,   41,   56], // 76%
        [  24,   33,   43,   57], // 77%
        [  26,   34,   45,   59], // 78%
        [  28,   37,   48,   61], // 79%
        [  30,   37,   48,   62], // 80%
        [  33,   42,   52,   64], // 81%
        [  35,   43,   53,   65], // 82%
        [  38,   47,   56,   68], // 83%
        [  40,   47,   56,   68], // 84%
        [  44,   53,   61,   71], // 85%
        [  45,   53,   62,   73], // 86%
        [  50,   58,   66,   75], // 87%
        [  51,   58,   66,   76], // 88%
        [  57,   63,   71,   79], // 89%
        [  60,   66,   72,   81], // 90%
        [  62,   68,   75,   83], // 91%
        [  70,   74,   80,   85], // 92%
        [  74,   78,   82,   88], // 93%
        [  85,   87,   90,   92], // 94%
        [  86,   88,   90,   92], // 95%
        [  87,   89,   91,   93], // 96%
        [  87,   89,   92,   94], // 97%
        [  89,   91,   93,   96], // 98%
        [  93,   94,   96,   97], // 99%
        [  94,   95,   97,   98], // 100%
    ],
    // cover 0.80
    [
        //  N=5   N=4   N=3   N=2
        [   0,    1,    4,   13], // 40%
        [   0,    1,    5,   13], // 41%
        [   0,    1,    5,   14], // 42%
        [   0,    2,    6,   15], // 43%
        [   0,    2,    6,   16], // 44%
        [   0,    2,    7,   17], // 45%
        [   1,    3,    8,   18], // 46%
        [   1,    4,    9,   20], // 47%
        [   1,    4,    9,   20], // 48%
        [   2,    5,   11,   22], // 49%
        [   2,    5,   11,   22], // 50%
        [   2,    6,   12,   24], // 51%
        [   3,    6,   13,   25], // 52%
        [   3,    7,   14,   26], // 53%
        [   4,    8,   14,   27], // 54%
        [   4,    8,   15,   28], // 55%
        [   5,    9,   17,   30], // 56%
        [   5,    9,   17,   30], // 57%
        [   6,   11,   19,   32], // 58%
        [   7,   12,   20,   34], // 59%
        [   8,   12,   20,   34], // 60%
        [   9,   14,   22,   37], // 61%
        [   9,   14,   23,   37], // 62%
        [  10,   16,   25,   39], // 63%
        [  11,   17,   26,   41], // 64%
        [  12,   18,   27,   41], // 65%
        [  13,   20,   28,   43], // 66%
        [  14,   21,   30,   45], // 67%
        [  15,   22,   31,   46], // 68%
        [  17,   24,   33,   48], // 69%
        [  17,   24,   34,   48], // 70%
        [  19,   26,   36,   50], // 71%
        [  20,   27,   37,   51], // 72%
        [  21,   29,   39,   53], // 73%
        [  23,   31,   41,   55], // 74%
        [  23,   31,   41,   55], // 75%
        [  26,   34,   44,   58], // 76%
        [  28,   36,   46,   59], // 77%
        [  29,   37,   47,   60], // 78%
        [  34,   41,   50,   62], // 79%
        [  34,   42,   51,   63], // 80%
        [  38,   45,   55,   66], // 81%
        [  39,   46,   55,   67], // 82%
        [  44,   51,   60,   70], // 83%
        [  44,   51,   60,   70], // 84%
        [  49,   56,   64,   73], // 85%
        [  50,   57,   64,   74], // 86%
        [  57,   63,   69,   77], // 87%
        [  58,   64,   70,   78], // 88%
        [  68,   71,   76,   82], // 89%
        [  68,   72,   77,   83], // 90%
        [  75,   79,   81,   85], // 91%
        [  86,   87,   89,   90], // 92%
        [  88,   89,   90,   92], // 93%
        [  90,   91,   92,   93], // 94%
        [  91,   92,   93,   94], // 95%
        [  92,   94,   94,   95], // 96%
        [  93,   94,   95,   96], // 97%
        [  94,   95,   95,   96], // 98%
        [  94,   95,   96,   98], // 99%
        [  95,   96,   97,   98], // 100%
    ],
    // cover 0.75
    [
        //  N=5   N=4   N=3   N=2
        [   1,    2,    6,   15], // 40%
        [   1,    3,    7,   16], // 41%
        [   1,    3,    8,   17], // 42%
        [   2,    4,    9,   18], // 43%
        [   2,    4,    9,   19], // 44%
        [   2,    5,   10,   20], // 45%
        [   3,    5,   10,   21], // 46%
        [   3,    6,   12,   22], // 47%
        [   3,    6,   12,   23], // 48%
        [   4,    8,   14,   25], // 49%
        [   4,    8,   14,   25], // 50%
        [   5,    8,   15,   26], // 51%
        [   5,    9,   16,   27], // 52%
        [   6,   10,   17,   29], // 53%
        [   6,   11,   18,   30], // 54%
        [   7,   11,   18,   31], // 55%
        [   8,   12,   20,   32], // 56%
        [   8,   13,   20,   33], // 57%
        [  10,   14,   22,   35], // 58%
        [  10,   15,   23,   37], // 59%
        [  11,   16,   24,   37], // 60%
        [  12,   18,   26,   39], // 61%
        [  13,   18,   26,   40], // 62%
        [  14,   20,   28,   42], // 63%
        [  16,   22,   30,   43], // 64%
        [  16,   22,   31,   44], // 65%
        [  17,   23,   32,   45], // 66%
        [  18,   25,   33,   47], // 67%
        [  19,   26,   35,   48], // 68%
        [  21,   27,   36,   50], // 69%
        [  22,   29,   37,   51], // 70%
        [  24,   30,   39,   52], // 71%
        [  25,   32,   41,   53], // 72%
        [  26,   33,   42,   55], // 73%
        [  29,   35,   44,   57], // 74%
        [  29,   36,   45,   57], // 75%
        [  32,   39,   48,   60], // 76%
        [  34,   41,   50,   61], // 77%
        [  36,   43,   51,   62], // 78%
        [  40,   46,   54,   65], // 79%
        [  40,   46,   54,   65], // 80%
        [  46,   52,   59,   68], // 81%
        [  46,   52,   60,   69], // 82%
        [  53,   59,   65,   73], // 83%
        [  54,   60,   66,   73], // 84%
        [  63,   67,   73,   78], // 85%
        [  68,   71,   75,   79], // 86%
        [  78,   80,   82,   85], // 87%
        [  79,   81,   83,   85], // 88%
        [  83,   85,   86,   87], // 89%
        [  85,   86,   87,   89], // 90%
        [  86,   88,   89,   90], // 91%
        [  88,   89,   90,   91], // 92%
        [  90,   90,   91,   92], // 93%
        [  91,   92,   92,   93], // 94%
        [  92,   93,   94,   94], // 95%
        [  94,   94,   95,   95], // 96%
        [  95,   95,   96,   96], // 97%
        [  95,   96,   97,   97], // 98%
        [  96,   96,   97,   98], // 99%
        [  97,   98,   98,   99], // 100%
    ],
];

/////////////////////////// END ALL ////////////////////////
