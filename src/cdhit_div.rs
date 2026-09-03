// =============================================================================
// CD-HIT-DIV
// http://cd-hit.org/
//
// Rust translation of `cdhit-div.c++`.
// =============================================================================

use crate::cdhit_common::*;
use crate::cdhit_utility::{current_time, print_usage_div};

////////////////////////////////////  MAIN /////////////////////////////////////
pub fn main(argv: &[String]) {
    let argc = argv.len();

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let alphabet = Alphabet::new();

    let mut db_in = String::new();
    let mut db_out = String::new();
    let mut div = 1;

    let begin_time = current_time();

    // ***********************************    parse command line and open file
    if argc < 5 {
        print_usage_div(&argv[0]);
    }
    let mut i = 1;
    while i < argc {
        match argv[i].as_str() {
            "-i" => {
                i += 1;
                db_in = argv[i].clone();
            }
            "-o" => {
                i += 1;
                db_out = argv[i].clone();
            }
            "-div" => {
                i += 1;
                div = argv[i].parse::<i32>().unwrap_or(0);
            }
            _ => print_usage_div(&argv[0]),
        }
        i += 1;
    }
    if div <= 1 {
        println!("Warning: -div must be greater than 1.");
        println!("Warning: no database is writen.");
        return;
    }
    options.store_disk = true;

    seq_db.read(&db_in, &options);
    println!("total seq: {}", seq_db.sequences.len());
    seq_db.sort_divide(&mut options, true, &alphabet);

    println!("writing new databases");
    seq_db.divide_save(&db_in, &db_out, div, &options);

    let end_time = current_time();
    println!("Total CPU time {:.2}", end_time - begin_time);

    println!("program completed !\n");
} // END int main
