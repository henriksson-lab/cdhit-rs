// =============================================================================
// CD-HIT-2D
// http://cd-hit.org/
//
// Rust translation of `cdhit-2d.c++`.
// =============================================================================

use cdhit_rs::cdhit_common::*;
use cdhit_rs::cdhit_utility::{current_time, print_usage_2d};

////////////////////////////////////  MAIN /////////////////////////////////////
fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let argc = argv.len();

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let mut seq_db2 = SequenceDB::default();
    let mut mat = ScoreMatrix::new();
    let mut alphabet = Alphabet::new();

    let begin_time = current_time();

    // ***********************************    parse command line and open file
    if argc < 7 {
        print_usage_2d(&argv[0]);
    }
    let num_procs = std::thread::available_parallelism()
        .map(|n| n.get() as i32)
        .unwrap_or(1);
    if !options.set_options(&argv, true, false, num_procs, &mut mat, &mut alphabet) {
        print_usage_2d(&argv[0]);
    }
    options.validate();

    let db_in = options.input.clone();
    let db_in2 = options.input2.clone();
    let db_out = options.output.clone();

    let naa_tab = Naa::new(MAX_UAA);
    options.naan = naa_tab.naan_array[options.naa as usize];
    seq_db.naan = naa_tab.naan_array[options.naa as usize];
    seq_db2.naan = naa_tab.naan_array[options.naa as usize];

    seq_db.read(&db_in, &options);
    println!("total seq in db1: {}", seq_db.sequences.len());

    seq_db2.read(&db_in2, &options);
    println!("total seq in db2: {}", seq_db2.sequences.len());

    seq_db.sort_divide(&mut options, true, &alphabet);
    seq_db2.sort_divide(&mut options, false, &alphabet);
    let opts = options.clone();
    seq_db2.cluster_to(&mut seq_db, &opts, &mat, &naa_tab, &[]);

    println!("writing non-redundant sequences from db2");
    seq_db2.write_clusters(&db_in2, &db_out, &options);

    // write a backup clstr file in case next step crashes
    seq_db2.write_extra_2d(&seq_db, &options);
    println!("program completed !\n");
    let end_time = current_time();
    println!("Total CPU time {:.2}", end_time - begin_time);
} // END int main
