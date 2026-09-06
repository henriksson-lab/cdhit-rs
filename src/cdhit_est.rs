// =============================================================================
// CD-HIT-EST
// http://cd-hit.org/
//
// Rust translation of `cdhit-est.c++`.
// =============================================================================

use crate::cdhit_common::*;
use crate::cdhit_utility::{current_time, print_usage_est};

/// `#undef MAX_UAA` / `#define MAX_UAA 4` -- the EST front-ends re-define the alphabet size
/// before calling `InitNAA`, so nucleotide words are base 4 rather than base 21.
const MAX_UAA_EST: i32 = 4;

////////////////////////////////////  MAIN /////////////////////////////////////
pub fn main(argv: &[String]) {
    let argc = argv.len();

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let mut mat = ScoreMatrix::new();
    let mut alphabet = Alphabet::new();

    options.cluster_thd = 0.95;
    options.naa = 10;
    options.naa_top_limit = 12;
    alphabet.setaa_to_na();
    mat.set_to_na(); //mat.set_gap(-6,-1);

    let begin_time = current_time();

    // ***********************************    parse command line and open file
    if argc < 5 {
        print_usage_est(&argv[0]);
    }
    let num_procs = std::thread::available_parallelism()
        .map(|n| n.get() as i32)
        .unwrap_or(1);
    if !options.set_options(argv, false, true, num_procs, &mut mat, &mut alphabet) {
        print_usage_est(&argv[0]);
    }
    options.validate();

    let db_in = options.input.clone();
    let db_in_pe = options.input_pe.clone();
    let db_out = options.output.clone();
    let db_out_pe = options.output_pe.clone();

    let naa_tab = Naa::new(MAX_UAA_EST);
    seq_db.naan = naa_tab.naan_array[options.naa as usize];

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

    println!("total seq: {}", seq_db.sequences.len());
    seq_db.sort_divide(&mut options, true, &alphabet);
    seq_db.do_clustering(&options, &mat, &naa_tab, &comp_aan_idx);

    println!("writing new database");
    if options.pe_mode != 0 {
        seq_db.write_clusters_pe(&db_in, &db_in_pe, &db_out, &db_out_pe, &options);
    } else {
        seq_db.write_clusters(&db_in, &db_out, &options);
    }

    // write a backup clstr file in case next step crashes
    seq_db.write_extra_1d(&options);
    println!("program completed !\n");
    let end_time = current_time();
    println!("Total CPU time {:.2}", end_time - begin_time);
}
