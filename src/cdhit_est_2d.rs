// =============================================================================
// CD-HIT-EST-2D
// http://cd-hit.org/
//
// Rust translation of `cdhit-est-2d.c++`.
// =============================================================================

use crate::cdhit_common::*;
use crate::cdhit_utility::{current_time, print_usage_est_2d};

/// `#undef MAX_UAA` / `#define MAX_UAA 4`: nucleotide words are base 4, not base 21.
const MAX_UAA_EST: i32 = 4;

////////////////////////////////////  MAIN /////////////////////////////////////
pub fn main(argv: &[String]) {
    let argc = argv.len();

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let mut seq_db2 = SequenceDB::default();
    let mut mat = ScoreMatrix::new();
    let mut alphabet = Alphabet::new();

    options.cluster_thd = 0.95;
    options.naa = 10;
    options.naa_top_limit = 12;
    alphabet.setaa_to_na();
    mat.set_to_na();

    let begin_time = current_time();

    // ***********************************    parse command line and open file
    if argc < 7 {
        print_usage_est_2d(&argv[0]);
    }
    let num_procs = std::thread::available_parallelism()
        .map(|n| n.get() as i32)
        .unwrap_or(1);
    if !options.set_options(argv, true, true, num_procs, &mut mat, &mut alphabet) {
        print_usage_est_2d(&argv[0]);
    }
    options.validate();

    let db_in = options.input.clone();
    let db_in_pe = options.input_pe.clone();
    let db_in2 = options.input2.clone();
    let db_in2_pe = options.input2_pe.clone();
    let db_out = options.output.clone();
    let db_out_pe = options.output_pe.clone();

    let naa_tab = Naa::new(MAX_UAA_EST);
    options.naan = naa_tab.naan_array[options.naa as usize];
    seq_db.naan = naa_tab.naan_array[options.naa as usize];
    seq_db2.naan = naa_tab.naan_array[options.naa as usize];

    let mut comp_aan_idx: Vec<i32> = Vec::new();
    if options.option_r != 0 {
        comp_aan_idx.resize(seq_db.naan as usize, 0);
        make_comp_short_word_index(options.naa, &naa_tab.naan_array, &mut comp_aan_idx);
    }

    if options.pe_mode != 0 {
        seq_db.read_pe(&db_in, &db_in_pe, &options);
    } else {
        seq_db.read(&db_in, &options);
    }
    println!("total seq in db1: {}", seq_db.sequences.len());

    if options.pe_mode != 0 {
        seq_db2.read_pe(&db_in2, &db_in2_pe, &options);
    } else {
        seq_db2.read(&db_in2, &options);
    }
    println!("total seq in db2: {}", seq_db2.sequences.len());

    seq_db.sort_divide(&mut options, true, &alphabet);
    seq_db2.sort_divide(&mut options, false, &alphabet);
    let opts = options.clone();
    seq_db2.cluster_to(&mut seq_db, &opts, &mat, &naa_tab, &comp_aan_idx);

    println!("writing non-redundant sequences from db2");
    if options.pe_mode != 0 {
        seq_db2.write_clusters_pe(&db_in2, &db_in2_pe, &db_out, &db_out_pe, &options);
    } else {
        seq_db2.write_clusters(&db_in2, &db_out, &options);
    }

    // write a backup clstr file in case next step crashes
    seq_db2.write_extra_2d(&seq_db, &options);
    println!("program completed !\n");
    let end_time = current_time();
    println!("Total CPU time {:.2}", end_time - begin_time);
} // END int main
