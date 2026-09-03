// =============================================================================
// cdhit -- the single front-end
//
// DEVIATION from CD-HIT 4.8.1, which ships six separate executables. They are subcommands
// here, all in one binary. The programs themselves are unchanged: each subcommand is the
// faithful translation of the corresponding `cdhit-*.c++` main, living in its own module.
//
// Two ways to select one, so existing scripts keep working:
//   * a subcommand      -- `cdhit est -i in.fna -o out`
//   * the program name  -- argv[0] matching a program, e.g. via a symlink named
//                          `cd-hit-est` pointing at this binary
// =============================================================================

use cdhit_rs::{program, PROGRAMS};

fn usage() -> ! {
    println!(
        "cdhit {} -- CD-HIT {}, translated to Rust\n",
        env!("CARGO_PKG_VERSION"),
        cdhit_rs::cdhit_common::CDHIT_VERSION
    );
    println!("Usage: cdhit <command> [options]\n");
    println!("Commands (upstream name, then the short form):\n");
    for (full, alias, _) in PROGRAMS {
        let what = match *full {
            "cd-hit" => "cluster a protein database",
            "cd-hit-est" => "cluster a nucleotide database, both strands",
            "cd-hit-2d" => "compare a protein db2 against db1, db1 unclustered",
            "cd-hit-est-2d" => "the same for nucleotides",
            "cd-hit-454" => "remove 454 duplicates, tolerating indels",
            "cd-hit-div" => "split a database into N parts",
            _ => "",
        };
        println!("   {:<14} {:<8} {}", full, alias, what);
    }
    println!("\n`cdhit <command> -h` prints that command's options.\n");
    println!("A symlink named after a command runs it directly, so scripts calling");
    println!("`cd-hit-est` keep working:  ln -s cdhit cd-hit-est\n");
    std::process::exit(1);
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();

    // Invoked through a name that is itself a program (a symlink, or a copy)?
    // `file_stem`, not `file_name`: on Windows a copy named `cd-hit-est.exe` must still
    // match, and for an extensionless name the two are the same.
    let arg0 = std::path::Path::new(&argv[0])
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    if let Some(run) = program(&arg0) {
        run(&argv);
        return;
    }

    if argv.len() < 2 {
        usage();
    }
    let Some(run) = program(&argv[1]) else {
        if matches!(argv[1].as_str(), "-h" | "--help" | "help") {
            usage();
        }
        eprintln!("cdhit: unknown command `{}`\n", argv[1]);
        usage();
    };

    // Hand the subcommand an argument vector shaped like the one the standalone program
    // would have seen: a program name followed by its own options. The name it prints in
    // its usage text is `cdhit <command>`, which is what the user would have to retype.
    let mut sub: Vec<String> = Vec::with_capacity(argv.len() - 1);
    sub.push(format!("{} {}", arg0, argv[1]));
    sub.extend_from_slice(&argv[2..]);
    run(&sub);
}
