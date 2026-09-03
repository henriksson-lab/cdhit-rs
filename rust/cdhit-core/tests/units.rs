//! Unit tests for the leaf functions of the translation -- the ones whose C++ behaviour is
//! subtle enough that a differential run over whole FASTA files would not pin it down.
//!
//! The end-to-end check is `tests/difftest.sh`, which diffs real output against the C++
//! binaries; these cover the corners that inputs rarely reach.

use cdhit_core::cdhit_common::*;

/// The C library's `atoi`: parse a leading integer, yield 0 when there isn't one.
#[test]
fn atoi_matches_c() {
    // `atoi` is private, so exercise it through the option parser, which is where it runs.
    let mut o = Options::default();
    let mut m = ScoreMatrix::new();
    let mut a = Alphabet::new();

    for (input, expect) in [
        ("12", 12),
        ("12abc", 12),
        ("abc", 0),
        ("", 0),
        ("-5", -5),
        ("+7", 7),
        ("  9", 9),
        ("3.7", 3),
    ] {
        assert!(o.set_option("-b", input, 1, &mut m, &mut a));
        assert_eq!(o.band_width, expect, "atoi({input:?})");
    }
}

/// The C library's `atof`.
#[test]
fn atof_matches_c() {
    let mut o = Options::default();
    let mut m = ScoreMatrix::new();
    let mut a = Alphabet::new();

    for (input, expect) in [
        ("0.9", 0.9),
        ("0.9xyz", 0.9),
        ("xyz", 0.0),
        ("1e-2", 0.01),
        (".5", 0.5),
        ("-0.25", -0.25),
    ] {
        assert!(o.set_option("-c", input, 1, &mut m, &mut a));
        assert!(
            (o.cluster_thd - expect).abs() < 1e-12,
            "atof({input:?}) = {} want {expect}",
            o.cluster_thd
        );
    }
}

/// `InitNAA` computes powers in `int`. For the protein alphabet everything from `NAA8` up
/// overflows; the EST front-ends avoid it by redefining `MAX_UAA` to 4 first. Values here
/// were taken from a C++ probe built against the vendored headers.
#[test]
fn naa_powers_wrap_like_int() {
    let protein = Naa::new(MAX_UAA);
    assert_eq!(protein.naan_array[1], 21);
    assert_eq!(protein.naan_array[5], 4084101);
    assert_eq!(protein.naan_array[4], 194481);
    assert_eq!(protein.naan_array[6], 85766121);
    assert_eq!(protein.naan_array[8], -831846303);
    assert_eq!(protein.naan_array[10], -1771999463);
    assert_eq!(protein.naan_array[12], 232284689);

    let dna = Naa::new(4);
    assert_eq!(dna.naan_array[4], 256);
    assert_eq!(dna.naan_array[8], 65536);
    assert_eq!(dna.naan_array[10], 1048576);
    assert_eq!(dna.naan_array[11], 4194304);
}

/// `Vector::Append` reserves `n + n/5 + 1`, and `MinimalMemory` reports the resulting
/// capacity, which feeds the table limits that decide where clustering blocks are cut.
#[test]
fn cpp_vector_growth_schedule() {
    // Hand-simulated from `Vector::Append`: reserve when `n + 1 >= capacity`.
    let expected = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 13, 15];
    for (n, want) in expected.iter().enumerate() {
        assert_eq!(cpp_vector_capacity(n), *want, "after {n} appends");
    }
    // Growth stays at roughly 1.2x, never doubling the way `Vec` does.
    let c = cpp_vector_capacity(1_000_000);
    assert!(
        (1_000_000..1_260_000).contains(&c),
        "capacity for 1e6 appends was {c}"
    );
}

/// `format_seq` upper-cases and drops everything that is not a letter.
#[test]
fn format_seq_keeps_letters_only() {
    let mut s = b"ac gt\n12-NNN*".to_vec();
    format_seq(&mut s);
    assert_eq!(s, b"ACGTNNN".to_vec());

    let mut empty = b"123 \n".to_vec();
    format_seq(&mut empty);
    assert!(empty.is_empty());
}

/// `Sequence::Format` rejects a record containing anything but letters and whitespace, and
/// strips one trailing `*` plus trailing whitespace.
#[test]
fn sequence_format_validates_and_strips() {
    let mut s = Sequence::default();
    s.assign(b"acgt\n");
    assert_eq!(s.format(), 0);
    assert_eq!(s.data, b"ACGT".to_vec());
    assert_eq!(s.size, 4);

    // a single trailing stop codon is dropped
    let mut stop = Sequence::default();
    stop.assign(b"MKV*\n");
    assert_eq!(stop.format(), 0);
    assert_eq!(stop.data, b"MKV".to_vec());

    // digits and punctuation make the record invalid; the count is returned
    let mut bad = Sequence::default();
    bad.assign(b"ACGT123\n");
    assert_eq!(bad.format(), 3);
}

/// `make_comp_iseq` reverse-complements in the encoded alphabet: A,C,G,T = 0..=3 swap in
/// pairs, N (4) and the mask code (5) are left alone.
#[test]
fn make_comp_iseq_reverse_complements() {
    let fwd = [0u8, 1, 2, 3, 4, 5]; // A C G T N X
    let mut rc = vec![0u8; fwd.len()];
    make_comp_iseq(fwd.len(), &mut rc, &fwd);
    assert_eq!(rc, vec![5, 4, 0, 1, 2, 3]);

    // applying it twice is the identity
    let mut back = vec![0u8; fwd.len()];
    make_comp_iseq(rc.len(), &mut back, &rc);
    assert_eq!(back, fwd.to_vec());
}

/// `make_comp_short_word_index` maps a word code to its reverse-complement word code.
#[test]
fn comp_word_index_is_an_involution() {
    let naa = 5;
    let naa_tab = Naa::new(4);
    let naan = naa_tab.naan_array[naa as usize] as usize;
    let mut idx = vec![0i32; naan];
    make_comp_short_word_index(naa, &naa_tab.naan_array, &mut idx);

    assert_eq!(idx.len(), 1024);
    for (code, &comp) in idx.iter().enumerate() {
        assert!((comp as usize) < naan, "code {code} maps out of range");
        assert_eq!(idx[comp as usize] as usize, code, "not an involution at {code}");
    }
    // AAAAA (all zeros) <-> TTTTT (all threes)
    assert_eq!(idx[0] as usize, naan - 1);
}

/// `quick_sort_idxr` sorts descending and carries a companion index array along. Cluster
/// output ordering under `-sc` / `-sf` depends on how it breaks ties, so pin that too.
#[test]
fn quick_sort_idxr_sorts_descending_with_companion() {
    let mut a = vec![3, 1, 4, 1, 5, 9, 2, 6];
    let mut idx: Vec<i32> = (0..a.len() as i32).collect();
    let n = a.len() as i32;
    quick_sort_idxr(&mut a, &mut idx, 0, n - 1);
    assert_eq!(a, vec![9, 6, 5, 4, 3, 2, 1, 1]);
    // the companion array must still point at the original positions
    let original = [3, 1, 4, 1, 5, 9, 2, 6];
    for (k, &i) in idx.iter().enumerate() {
        assert_eq!(original[i as usize], a[k]);
    }

    let mut single = vec![7];
    let mut si = vec![0];
    quick_sort_idxr(&mut single, &mut si, 0, 0);
    assert_eq!(single, vec![7]);
}

/// `quick_sort_idx` is the ascending counterpart.
#[test]
fn quick_sort_idx_sorts_ascending() {
    let mut a = vec![3, 1, 4, 1, 5, 9, 2, 6];
    let mut idx: Vec<i32> = (0..a.len() as i32).collect();
    let n = a.len() as i32;
    quick_sort_idx(&mut a, &mut idx, 0, n - 1);
    assert_eq!(a, vec![1, 1, 2, 3, 4, 5, 6, 9]);
}

/// Every score is scaled by `MAX_SEQ`, which is what makes `local_band_align`'s
/// `extra_score` bonus a tie-break rather than a contribution -- and what forces the DP
/// cells to be 64-bit.
#[test]
fn score_matrix_is_scaled_by_max_seq() {
    let m = ScoreMatrix::new();
    assert_eq!(m.matrix[0][0], 4 * MAX_SEQ); // BLOSUM62 A/A
    assert_eq!(m.matrix[0][1], -MAX_SEQ); // A/R
    assert_eq!(m.matrix[1][0], -MAX_SEQ); // mirrored
    assert_eq!(m.gap, -11 * MAX_SEQ);
    assert_eq!(m.ext_gap, -MAX_SEQ);

    let mut na = ScoreMatrix::new();
    na.set_to_na();
    assert_eq!(na.matrix[0][0], 2 * MAX_SEQ); // A/A
    assert_eq!(na.matrix[0][1], -2 * MAX_SEQ); // A/C
    assert_eq!(na.matrix[3][4], MAX_SEQ); // T/U scores as a match
    assert_eq!(na.gap, -6 * MAX_SEQ);

    // `-mismatch` rewrites every off-diagonal cell but restores the T/U match.
    na.set_mismatch(-3);
    assert_eq!(na.matrix[0][1], -3 * MAX_SEQ);
    assert_eq!(na.matrix[3][4], MAX_SEQ);
    assert_eq!(na.matrix[0][0], 2 * MAX_SEQ); // diagonal untouched
}

/// `setaa_to_na` replaces the protein alphabet with the nucleotide one in place.
#[test]
fn setaa_to_na_switches_alphabet() {
    let mut a = Alphabet::new();
    assert_eq!(a.aa2idx[(b'M' - b'A') as usize], 12); // methionine
    a.setaa_to_na();
    for i in 0..26 {
        assert_eq!(a.aa2idx[i], a.na2idx[i]);
    }
    assert_eq!(a.aa2idx[(b'A' - b'A') as usize], 0);
    assert_eq!(a.aa2idx[(b'C' - b'A') as usize], 1);
    assert_eq!(a.aa2idx[(b'G' - b'A') as usize], 2);
    assert_eq!(a.aa2idx[(b'T' - b'A') as usize], 3);
    assert_eq!(a.aa2idx[(b'U' - b'A') as usize], 3);
    assert_eq!(a.aa2idx[(b'N' - b'A') as usize], 4);
}

/// `upper_bound_length_rep` turns the length-difference options into the cap that lets the
/// length-sorted scan stop early.
#[test]
fn upper_bound_length_rep_combines_options() {
    let mut o = Options::default();
    // no constraints at all: effectively unbounded
    assert_eq!(upper_bound_length_rep(100, &o), 99999999);

    o.diff_cutoff = 0.8; // -s
    assert_eq!(upper_bound_length_rep(100, &o), 125);

    o.diff_cutoff_aa = 10; // -S, tighter here
    assert_eq!(upper_bound_length_rep(100, &o), 110);

    o.long_coverage = 0.9; // -aL beats -s
    assert_eq!(upper_bound_length_rep(100, &o), 110);
    o.diff_cutoff_aa = 99999999;
    assert_eq!(upper_bound_length_rep(100, &o), 111);
}

/// `cal_aax_cutoff` derives the word-count thresholds; with `-t 0` it uses the closed form,
/// otherwise the empirical `naa_stat` table can only raise them.
#[test]
fn cal_aax_cutoff_tolerance_zero_is_closed_form() {
    let (mut a1, mut a2, mut an) = (0.0, 0.0, 0.0);
    cal_aax_cutoff(&mut a1, &mut a2, &mut an, 0.9, 0, NAA_STAT_START_PERCENT, &NAA_STAT, 5);
    assert!((a1 - 0.9).abs() < 1e-12);
    assert!((a2 - 0.8).abs() < 1e-12); // 1 - (1-0.9)*2
    assert!((an - 0.5).abs() < 1e-12); // 1 - (1-0.9)*5

    let (mut b1, mut b2, mut bn) = (0.0, 0.0, 0.0);
    cal_aax_cutoff(&mut b1, &mut b2, &mut bn, 0.9, 2, NAA_STAT_START_PERCENT, &NAA_STAT, 5);
    assert!((b1 - a1).abs() < 1e-12, "aa1 is never adjusted");
    assert!(b2 >= a2, "the table can only raise aa2");
    assert!(bn >= an, "the table can only raise aan");
}

/// The empirical table is 5 tolerances x 61 identity percentages x 4 word lengths, and is
/// indexed as `[tolerance-1][pct-40][5-word_len]`.
#[test]
fn naa_stat_table_shape() {
    assert_eq!(NAA_STAT.len(), 5);
    assert_eq!(NAA_STAT[0].len(), 61);
    assert_eq!(NAA_STAT[0][0].len(), 4);
    assert_eq!(NAA_STAT_START_PERCENT, 40);
    // The 100%-identity row of each tolerance block, straight from `cdhit-common.c++`.
    assert_eq!(NAA_STAT[0][60], [91, 93, 95, 97]); // cover 0.99
    assert_eq!(NAA_STAT[1][60], [92, 93, 95, 97]); // cover 0.95
    assert_eq!(NAA_STAT[2][60], [94, 95, 97, 98]); // cover 0.90
    assert_eq!(NAA_STAT[3][60], [95, 96, 97, 98]); // cover 0.80
    assert_eq!(NAA_STAT[4][60], [97, 98, 98, 99]); // cover 0.75
    // and the 40% row, where a 5-mer requirement has not yet appeared
    assert_eq!(NAA_STAT[0][0], [0, 0, 0, 7]);
    // within a tolerance, the requirement rises with identity
    for (t, block) in NAA_STAT.iter().enumerate() {
        for (w, (&hi, &lo)) in block[60].iter().zip(block[0].iter()).enumerate() {
            assert!(hi >= lo, "tolerance {t}, column {w} is not monotone");
        }
    }
}

/// `Sequence::clear` drops the residues but keeps `size`, which the writer still needs.
#[test]
fn sequence_clear_keeps_size() {
    let mut s = Sequence::default();
    s.assign(b"ACGTACGT");
    assert_eq!(s.size, 8);
    s.clear();
    assert_eq!(s.size, 8, "size must survive for the .clstr output");
    assert!(s.data.is_empty());
}

/// `Sequence::trim` implements `-cx` / `-cy`, and is a no-op when the sequence is already
/// short enough.
#[test]
fn sequence_trim_is_a_prefix() {
    let mut s = Sequence::default();
    s.assign(b"ACGTACGTAC");
    s.trim(4);
    assert_eq!(s.data, b"ACGT".to_vec());
    assert_eq!(s.size, 4);
    s.trim(100);
    assert_eq!(s.size, 4);
}

/// `Sequence::ConvertBases` rewrites letters to alphabet indices in place.
#[test]
fn convert_bases_encodes_in_place() {
    let mut a = Alphabet::new();
    a.setaa_to_na();
    let mut s = Sequence::default();
    s.assign(b"ACGTN");
    s.convert_bases(&a);
    assert_eq!(s.data, vec![0, 1, 2, 3, 4]);
}
