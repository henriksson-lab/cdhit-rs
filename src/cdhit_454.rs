// =============================================================================
// CD-HIT-454
// http://cd-hit.org/
//
// Rust translation of `cdhit-454.c++`.
// =============================================================================

use crate::cdhit_common::*;
use crate::cdhit_utility::{current_time, print_usage_454};

/// `#undef MAX_UAA` / `#define MAX_UAA 4`: nucleotide words are base 4, not base 21.
const MAX_UAA_EST: i32 = 4;

// matrix below is for highly similar seqs
// with this matrix, alignment can be maintained if there is a mismatch
// at very end of alignment
// |x||||||||
//  ^
//  mismatch will cost -2 with BLOSUM62_na, therefore the first 2 bases won't be in alignment
//  but with BLOSUM62_na2, the first 2 bases are in alignment
/// `myBLOSUM62_na2`
static MY_BLOSUM62_NA2: [i32; 21] = [
    2, // A
    -1, 2, // C
    -1, -1, 2, // G
    -1, -1, -1, 2, // T
    -1, -1, -1, 2, 2, // U
    -1, -1, -1, -1, -1, 2, // N
       //A  C  G  T  U  N
       //0  1  2  3  3  4
];

////////////////////////////////////  MAIN /////////////////////////////////////
pub fn main(argv: &[String]) {
    let argc = argv.len();

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let mut mat = ScoreMatrix::new();
    let mut alphabet = Alphabet::new();

    options.is_est = true;
    options.is_454 = true;
    options.naa = 10;
    options.naa_top_limit = 12;

    options.cluster_thd = 0.98;
    options.band_width = 10;
    options.print = 1;
    options.des_len = 0;
    options.option_r = 0;

    alphabet.setaa_to_na();
    mat.set_gap(-3, -1); //instead of -6 -1 to maintain maxium length of alignment
    mat.set_matrix(&MY_BLOSUM62_NA2);

    let begin_time = current_time();

    // ***********************************    parse command line and open file
    if argc < 5 {
        print_usage_454(&argv[0]);
    }
    let num_procs = std::thread::available_parallelism()
        .map(|n| n.get() as i32)
        .unwrap_or(1);
    if !options.set_options(argv, false, true, num_procs, &mut mat, &mut alphabet) {
        print_usage_454(&argv[0]);
    }
    options.validate();

    let db_in = options.input.clone();
    let db_out = options.output.clone();

    let naa_tab = Naa::new(MAX_UAA_EST);
    options.naan = naa_tab.naan_array[options.naa as usize];
    seq_db.naan = naa_tab.naan_array[options.naa as usize];

    seq_db.read(&db_in, &options);
    println!("total seq: {}", seq_db.sequences.len());
    seq_db.sort_divide(&mut options, true, &alphabet);

    seq_db.do_clustering(&options, &mat, &naa_tab, &[]);

    println!("writing new database");
    seq_db.write_clusters(&db_in, &db_out, &options);

    // write a backup clstr file in case next step crashes
    seq_db.write_extra_1d(&options);
    println!("program completed !\n");
    let end_time = current_time();
    println!("Total CPU time {:.2}", end_time - begin_time);
} // END int main
