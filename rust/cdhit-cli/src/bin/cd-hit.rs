// =============================================================================
// CD-HIT
// http://cd-hit.org/
//
// Rust translation of `cdhit.c++`.
// =============================================================================

use cdhit_core::cdhit_common::*;
use cdhit_core::cdhit_utility::{current_time, print_usage};

////////////////////////////////////  MAIN /////////////////////////////////////
fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let argc = argv.len();

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let mut mat = ScoreMatrix::new();
    let mut alphabet = Alphabet::new();

    let begin_time = current_time();

    // ***********************************    parse command line and open file
    if argc < 5 {
        print_usage(&argv[0]);
    }
    let num_procs = std::thread::available_parallelism()
        .map(|n| n.get() as i32)
        .unwrap_or(1);
    if !options.set_options(&argv, false, false, num_procs, &mut mat, &mut alphabet) {
        print_usage(&argv[0]);
    }
    options.validate();

    let db_in = options.input.clone();
    let db_out = options.output.clone();

    let naa_tab = Naa::new(MAX_UAA);
    options.naan = naa_tab.naan_array[options.naa as usize];
    seq_db.naan = naa_tab.naan_array[options.naa as usize];

    seq_db.read(&db_in, &options);
    println!("total seq: {}", seq_db.sequences.len());
    seq_db.sort_divide(&mut options, true, &alphabet);

    seq_db.do_clustering(&options.clone(), &mat, &naa_tab, &[]);

    println!("writing new database");
    seq_db.write_clusters(&db_in, &db_out, &options);

    // write a backup clstr file in case next step crashes
    seq_db.write_extra_1d(&options);
    println!("program completed !\n");
    let end_time = current_time();
    println!("Total CPU time {:.2}", end_time - begin_time);
} // END int main
