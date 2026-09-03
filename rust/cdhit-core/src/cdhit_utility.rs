//! Rust translation of `cdhit-utility.h` + `cdhit-utility.c++`.
//!
//! The usage text is transcribed mechanically from the C++ string literals, C line splices
//! resolved, so the help output is byte-identical apart from the build date.
//!
//! The usage-text constants keep their C++ identifiers verbatim: `-g` and `-G` (and `-s` /
//! `-S`, `-t` / `-T`, ...) are distinct options differing only in letter case, which
//! SCREAMING_SNAKE_CASE would collapse into one name.

#![allow(non_upper_case_globals)]

use crate::cdhit_common::{build_date, CDHIT_VERSION};

/// `current_time`: `((float)clock())/CLOCKS_PER_SEC`, i.e. CPU time consumed by the
/// process so far, summed over threads.
///
/// `libc` does not expose `clock()`, so this reads `CLOCK_PROCESS_CPUTIME_ID` directly;
/// that is the clock `clock()` is specified to report.
pub fn current_time() -> f32 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // SAFETY: `ts` is a valid, fully initialised `timespec` for the duration of the call.
    let rc = unsafe { libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut ts) };
    if rc != 0 {
        return 0.0;
    }
    ts.tv_sec as f32 + ts.tv_nsec as f32 / 1e9
}

// information
/// `cd_hit_ver`
pub fn cd_hit_ver() -> String {
    format!(
        "\t\t====== CD-HIT version {} (built on {}) ======",
        CDHIT_VERSION,
        build_date()
    )
}

/// `cd_hit_ref1`
pub const cd_hit_ref1: &str = "\"CD-HIT: a fast program for clustering and comparing large sets of protein or nucleotide sequences\", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659";

/// `cd_hit_ref2`
pub const cd_hit_ref2: &str = "\"CD-HIT: accelerated for clustering the next generation sequencing data\", Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. Bioinformatics, (2012) 28:3150-3152";

/// `cd_hit_ref3`
pub const cd_hit_ref3: &str = "\"Beifang Niu, Limin Fu, Shulei Sun and Weizhong Li. Artificial and natural duplicates in pyrosequencing reads of metagenomic data. BMC Bioinformatics (2010) 11:187";

/// `contacts`
pub const contacts: &str = "   Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n   For updated versions and information, please visit: http://cd-hit.org\n                                                    or https://github.com/weizhongli/cdhit\n\n   cd-hit web server is also available from http://cd-hit.org\n\n   If you find cd-hit useful, please kindly cite:\n\n";

/// `txt_option_i`
pub const txt_option_i: &str = "\tinput filename in fasta format, required, can be in .gz format\n";

/// `txt_option_j`
pub const txt_option_j: &str = "\tinput filename in fasta/fastq format for R2 reads if input are paired end (PE) files\n \t -i R1.fq -j R2.fq -o output_R1 -op output_R2 or\n \t -i R1.fa -j R2.fa -o output_R1 -op output_R2 \n";

/// `txt_option_i_2d`
pub const txt_option_i_2d: &str = "\tinput filename for db1 in fasta format, required, can be in .gz format\n";

/// `txt_option_i2`
pub const txt_option_i2: &str = "\tinput filename for db2 in fasta format, required, can be in .gz format\n";

/// `txt_option_j2`
pub const txt_option_j2: &str = "\tinput filename in fasta/fastq format for R2 reads if input are paired end (PE) files\n \t -i db1-R1.fq -j db1-R2.fq -i2 db2-R1.fq -j2 db2-R2.fq -o output_R1 -op output_R2 or\n \t -i db1-R1.fa -j db1-R2.fa -i2 db2-R1.fq -j2 db2-R2.fq -o output_R1 -op output_R2 \n";

/// `txt_option_o`
pub const txt_option_o: &str = "\toutput filename, required\n";

/// `txt_option_op`
pub const txt_option_op: &str = "\toutput filename for R2 reads if input are paired end (PE) files\n";

/// `txt_option_c`
pub const txt_option_c: &str = "\tsequence identity threshold, default 0.9\n \tthis is the default cd-hit's \"global sequence identity\" calculated as:\n \tnumber of identical amino acids or bases in alignment\n \tdivided by the full length of the shorter sequence\n";

/// `txt_option_G`
pub const txt_option_G: &str = "\tuse global sequence identity, default 1\n \tif set to 0, then use local sequence identity, calculated as :\n \tnumber of identical amino acids or bases in alignment\n \tdivided by the length of the alignment\n \tNOTE!!! don't use -G 0 unless you use alignment coverage controls\n \tsee options -aL, -AL, -aS, -AS\n";

/// `txt_option_g`
pub const txt_option_g: &str = "\t1 or 0, default 0\n \tby cd-hit's default algorithm, a sequence is clustered to the first \n \tcluster that meet the threshold (fast cluster). If set to 1, the program\n \twill cluster it into the most similar cluster that meet the threshold\n \t(accurate but slow mode)\n \tbut either 1 or 0 won't change the representatives of final clusters\n";

/// `txt_option_b`
pub const txt_option_b: &str = "\tband_width of alignment, default 20\n";

/// `txt_option_M`
pub const txt_option_M: &str = "\tmemory limit (in MB) for the program, default 800; 0 for unlimitted;\n";

/// `txt_option_n`
pub const txt_option_n: &str = "\tword_length, default 5, see user's guide for choosing it\n";

/// `txt_option_n_est`
pub const txt_option_n_est: &str = "\tword_length, default 10, see user's guide for choosing it\n";

/// `txt_option_l`
pub const txt_option_l: &str = "\tlength of throw_away_sequences, default 10\n";

/// `txt_option_t`
pub const txt_option_t: &str = "\ttolerance for redundance, default 2\n";

/// `txt_option_T`
pub const txt_option_T: &str = "\tnumber of threads, default 1; with 0, all CPUs will be used\n";

/// `txt_option_d`
pub const txt_option_d: &str = "\tlength of description in .clstr file, default 20\n \tif set to 0, it takes the fasta defline and stops at first space\n";

/// `txt_option_s`
pub const txt_option_s: &str = "\tlength difference cutoff, default 0.0\n \tif set to 0.9, the shorter sequences need to be\n \tat least 90% length of the representative of the cluster\n";

/// `txt_option_S`
pub const txt_option_S: &str = "\tlength difference cutoff in amino acid, default 999999\n \tif set to 60, the length difference between the shorter sequences\n \tand the representative of the cluster can not be bigger than 60\n";

/// `txt_option_s2`
pub const txt_option_s2: &str = "\tlength difference cutoff for db1, default 1.0\n \tby default, seqs in db1 >= seqs in db2 in a same cluster\n \tif set to 0.9, seqs in db1 may just >= 90% seqs in db2\n";

/// `txt_option_S2`
pub const txt_option_S2: &str = "\tlength difference cutoff, default 0\n \tby default, seqs in db1 >= seqs in db2 in a same cluster\n \tif set to 60, seqs in db2 may 60aa longer than seqs in db1\n";

/// `txt_option_aL`
pub const txt_option_aL: &str = "\talignment coverage for the longer sequence, default 0.0\n \tif set to 0.9, the alignment must covers 90% of the sequence\n";

/// `txt_option_AL`
pub const txt_option_AL: &str = "\talignment coverage control for the longer sequence, default 99999999\n \tif set to 60, and the length of the sequence is 400,\n \tthen the alignment must be >= 340 (400-60) residues\n";

/// `txt_option_aS`
pub const txt_option_aS: &str = "\talignment coverage for the shorter sequence, default 0.0\n \tif set to 0.9, the alignment must covers 90% of the sequence\n";

/// `txt_option_AS`
pub const txt_option_AS: &str = "\talignment coverage control for the shorter sequence, default 99999999\n \tif set to 60, and the length of the sequence is 400,\n \tthen the alignment must be >= 340 (400-60) residues\n";

/// `txt_option_A`
pub const txt_option_A: &str = "\tminimal alignment coverage control for the both sequences, default 0\n \talignment must cover >= this value for both sequences \n";

/// `txt_option_B`
pub const txt_option_B: &str = "\t1 or 0, default 0, by default, sequences are stored in RAM\n \tif set to 1, sequence are stored on hard drive\n \t!! No longer supported !!\n";

/// `txt_option_P`
pub const txt_option_P: &str = "\tinput paired end (PE) reads, default 0, single file\n \tif set to 1, please use -i R1 -j R2 to input both PE files\n";

/// `txt_option_cx`
pub const txt_option_cx: &str = "\tlength to keep after trimming the tail of sequence, default 0, not trimming\n \tif set to 50, the program only uses the first 50 letters of input sequence\n";

/// `txt_option_cy`
pub const txt_option_cy: &str = "\tlength to keep after trimming the tail of R2 sequence, default 0, not trimming\n \tif set to 50, the program only uses the first 50 letters of input R2 sequence\n \te.g. -cx 100 -cy 80 for paired end reads\n";

/// `txt_option_ap`
pub const txt_option_ap: &str = "\talignment position constrains,  default 0, no constrain\n \tif set to 1, the program will force sequences to align at beginings\n \twhen set to 1, the program only does +/+ alignment\n";

/// `txt_option_uL`
pub const txt_option_uL: &str = "\tmaximum unmatched percentage for the longer sequence, default 1.0\n \tif set to 0.1, the unmatched region (excluding leading and tailing gaps)\n \tmust not be more than 10% of the sequence\n";

/// `txt_option_uS`
pub const txt_option_uS: &str = "\tmaximum unmatched percentage for the shorter sequence, default 1.0\n \tif set to 0.1, the unmatched region (excluding leading and tailing gaps)\n \tmust not be more than 10% of the sequence\n";

/// `txt_option_U`
pub const txt_option_U: &str = "\tmaximum unmatched length, default 99999999\n \tif set to 10, the unmatched region (excluding leading and tailing gaps)\n \tmust not be more than 10 bases\n";

/// `txt_option_p`
pub const txt_option_p: &str = "\t1 or 0, default 0\n \tif set to 1, print alignment overlap in .clstr file\n";

/// `txt_option_r`
pub const txt_option_r: &str = "\t1 or 0, default 1, by default do both +/+ & +/- alignments\n \tif set to 0, only +/+ strand alignment\n";

/// `txt_option_bak`
pub const txt_option_bak: &str = "\twrite backup cluster file (1 or 0, default 0)\n";

/// `txt_option_sc`
pub const txt_option_sc: &str = "\tsort clusters by size (number of sequences), default 0, output clusters by decreasing length\n \tif set to 1, output clusters by decreasing size\n";

/// `txt_option_sf`
pub const txt_option_sf: &str = "\tsort fasta/fastq by cluster size (number of sequences), default 0, no sorting\n \tif set to 1, output sequences by decreasing cluster size\n \tthis can be very slow if the input is in .gz format\n";

/// `txt_option_mask`
pub const txt_option_mask: &str = "\tmasking letters (e.g. -mask NX, to mask out both 'N' and 'X')\n";

/// `txt_option_match`
pub const txt_option_match: &str = "\tmatching score, default 2 (1 for T-U and N-N)\n";

/// `txt_option_match2`
pub const txt_option_match2: &str = "\tmatching score, default 2\n";

/// `txt_option_mismatch`
pub const txt_option_mismatch: &str = "\tmismatching score, default -2\n";

/// `txt_option_mismatch2`
pub const txt_option_mismatch2: &str = "\tmismatching score, default -1\n";

/// `txt_option_gap`
pub const txt_option_gap: &str = "\tgap opening score, default -6\n";

/// `txt_option_gap2`
pub const txt_option_gap2: &str = "\tgap opening score, default -3\n";

/// `txt_option_gap_ext`
pub const txt_option_gap_ext: &str = "\tgap extension score, default -1\n";

/// `mytxt_option_c`
pub const mytxt_option_c: &str = "\tsequence identity threshold, default 0.98\n \tthis is a \"global sequence identity\" calculated as :\n \tnumber of identical amino acids or bases in alignment\n \tdivided by the full length of the shorter sequence + gaps\n";

/// `mytxt_option_b`
pub const mytxt_option_b: &str = "\tband_width of alignment, default 10\n";

/// `mytxt_option_n_est`
pub const mytxt_option_n_est: &str = "\tword_length, default 10, see user's guide for choosing it\n";

/// `mytxt_option_D`
pub const mytxt_option_D: &str = "\tmax size per indel, default 1\n";

/// `print_usage`
pub fn print_usage(arg: &str) -> ! {
    println!("{}\n", cd_hit_ver());
    println!("Usage: {} [Options] \n\nOptions\n", arg);
    print!("   -i{}", txt_option_i);
    print!("   -o{}", txt_option_o);
    print!("   -c{}", txt_option_c);
    print!("   -G{}", txt_option_G);
    print!("   -b{}", txt_option_b);
    print!("   -M{}", txt_option_M);
    print!("   -T{}", txt_option_T);
    print!("   -n{}", txt_option_n);
    print!("   -l{}", txt_option_l);
    print!("   -t{}", txt_option_t);
    print!("   -d{}", txt_option_d);
    print!("   -s{}", txt_option_s);
    print!("   -S{}", txt_option_S);
    print!("   -aL{}", txt_option_aL);
    print!("   -AL{}", txt_option_AL);
    print!("   -aS{}", txt_option_aS);
    print!("   -AS{}", txt_option_AS);
    print!("   -A{}", txt_option_A);
    print!("   -uL{}", txt_option_uL);
    print!("   -uS{}", txt_option_uS);
    print!("   -U{}", txt_option_U);
    print!("   -B{}", txt_option_B);
    print!("   -p{}", txt_option_p);
    print!("   -g{}", txt_option_g);
    print!("   -sc{}", txt_option_sc);
    print!("   -sf{}", txt_option_sf);
    print!("   -bak{}", txt_option_bak);
    println!("   -h\tprint this help\n");
    print!("{}", contacts);
    println!("   {}", cd_hit_ref1);
    println!("   {}\n\n", cd_hit_ref2);
    std::process::exit(1);
} // END print_usage

/// `print_usage_2d`
pub fn print_usage_2d(arg: &str) -> ! {
    println!("{}\n", cd_hit_ver());
    println!("Usage: {} [Options] \n\nOptions\n", arg);
    print!("   -i{}", txt_option_i_2d);
    print!("   -i2{}", txt_option_i2);
    print!("   -o{}", txt_option_o);
    print!("   -c{}", txt_option_c);
    print!("   -G{}", txt_option_G);
    print!("   -b{}", txt_option_b);
    print!("   -M{}", txt_option_M);
    print!("   -T{}", txt_option_T);
    print!("   -n{}", txt_option_n);
    print!("   -l{}", txt_option_l);
    print!("   -t{}", txt_option_t);
    print!("   -d{}", txt_option_d);
    print!("   -s{}", txt_option_s);
    print!("   -S{}", txt_option_S);
    print!("   -s2{}", txt_option_s2);
    print!("   -S2{}", txt_option_S2);
    print!("   -aL{}", txt_option_aL);
    print!("   -AL{}", txt_option_AL);
    print!("   -aS{}", txt_option_aS);
    print!("   -AS{}", txt_option_AS);
    print!("   -A{}", txt_option_A);
    print!("   -uL{}", txt_option_uL);
    print!("   -uS{}", txt_option_uS);
    print!("   -U{}", txt_option_U);
    print!("   -B{}", txt_option_B);
    print!("   -p{}", txt_option_p);
    print!("   -g{}", txt_option_g);
    print!("   -bak{}", txt_option_bak);
    println!("   -h\tprint this help\n");
    println!("   Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n");
    println!("   If you find cd-hit useful, please kindly cite:\n");
    println!("   {}", cd_hit_ref1);
    println!("   {}\n\n", cd_hit_ref2);
    std::process::exit(1);
} // END print_usage_2d

/// `print_usage_est`
pub fn print_usage_est(arg: &str) -> ! {
    println!("{}\n", cd_hit_ver());
    println!("Usage: {} [Options] \n\nOptions\n", arg);
    print!("   -i{}", txt_option_i);
    print!("   -j{}", txt_option_j);
    print!("   -o{}", txt_option_o);
    print!("   -op{}", txt_option_op);
    print!("   -c{}", txt_option_c);
    print!("   -G{}", txt_option_G);
    print!("   -b{}", txt_option_b);
    print!("   -M{}", txt_option_M);
    print!("   -T{}", txt_option_T);
    print!("   -n{}", txt_option_n_est);
    print!("   -l{}", txt_option_l);
    print!("   -d{}", txt_option_d);
    print!("   -s{}", txt_option_s);
    print!("   -S{}", txt_option_S);
    print!("   -aL{}", txt_option_aL);
    print!("   -AL{}", txt_option_AL);
    print!("   -aS{}", txt_option_aS);
    print!("   -AS{}", txt_option_AS);
    print!("   -A{}", txt_option_A);
    print!("   -uL{}", txt_option_uL);
    print!("   -uS{}", txt_option_uS);
    print!("   -U{}", txt_option_U);
    print!("   -B{}", txt_option_B);
    print!("   -P{}", txt_option_P);
    print!("   -cx{}", txt_option_cx);
    print!("   -cy{}", txt_option_cy);
    print!("   -ap{}", txt_option_ap);
    print!("   -p{}", txt_option_p);
    print!("   -g{}", txt_option_g);
    print!("   -r{}", txt_option_r);
    print!("   -mask{}", txt_option_mask);
    print!("   -match{}", txt_option_match);
    print!("   -mismatch{}", txt_option_mismatch);
    print!("   -gap{}", txt_option_gap);
    print!("   -gap-ext{}", txt_option_gap_ext);
    print!("   -bak{}", txt_option_bak);
    print!("   -sc{}", txt_option_sc);
    print!("   -sf{}", txt_option_sf);
    println!("   -h\tprint this help\n");
    print!("{}", contacts);
    println!("   {}", cd_hit_ref1);
    println!("   {}\n\n", cd_hit_ref2);
    std::process::exit(1);
} // END print_usage_est

/// `print_usage_est_2d`
pub fn print_usage_est_2d(arg: &str) -> ! {
    println!("{}\n", cd_hit_ver());
    println!("Usage: {} [Options] \n\nOptions\n", arg);
    print!("   -i{}", txt_option_i_2d);
    print!("   -i2{}", txt_option_i2);
    print!("   -j, -j2{}", txt_option_j2);
    print!("   -o{}", txt_option_o);
    print!("   -op{}", txt_option_op);
    print!("   -c{}", txt_option_c);
    print!("   -G{}", txt_option_G);
    print!("   -b{}", txt_option_b);
    print!("   -M{}", txt_option_M);
    print!("   -T{}", txt_option_T);
    print!("   -n{}", txt_option_n_est);
    print!("   -l{}", txt_option_l);
    print!("   -d{}", txt_option_d);
    print!("   -s{}", txt_option_s);
    print!("   -S{}", txt_option_S);
    print!("   -s2{}", txt_option_s2);
    print!("   -S2{}", txt_option_S2);
    print!("   -aL{}", txt_option_aL);
    print!("   -AL{}", txt_option_AL);
    print!("   -aS{}", txt_option_aS);
    print!("   -AS{}", txt_option_AS);
    print!("   -A{}", txt_option_A);
    print!("   -uL{}", txt_option_uL);
    print!("   -uS{}", txt_option_uS);
    print!("   -U{}", txt_option_U);
    print!("   -B{}", txt_option_B);
    print!("   -P{}", txt_option_P);
    print!("   -cx{}", txt_option_cx);
    print!("   -cy{}", txt_option_cy);
    print!("   -p{}", txt_option_p);
    print!("   -g{}", txt_option_g);
    print!("   -r{}", txt_option_r);
    print!("   -mask{}", txt_option_mask);
    print!("   -match{}", txt_option_match);
    print!("   -mismatch{}", txt_option_mismatch);
    print!("   -gap{}", txt_option_gap);
    print!("   -gap-ext{}", txt_option_gap_ext);
    print!("   -bak{}", txt_option_bak);
    println!("   -h\tprint this help\n");
    print!("{}", contacts);
    println!("   {}", cd_hit_ref1);
    println!("   {}\n\n", cd_hit_ref2);
    std::process::exit(1);
} // END print_usage_est_2d

/// `print_usage_div`
pub fn print_usage_div(arg: &str) -> ! {
    println!("{}\n", cd_hit_ver());
    println!("Usage: {} [Options] \n\nOptions\n", arg);
    println!("Options \n");
    println!("   -i in_dbname, required");
    println!("   -o out_dbname, required");
    println!("   -div number of divide, required ");
    std::process::exit(1);
} // END print_usage_div

/// `print_usage_454`
pub fn print_usage_454(arg: &str) -> ! {
    println!("{}\n", cd_hit_ver());
    println!("Usage: {} [Options] \n\nOptions\n", arg);
    print!("   -i{}", txt_option_i);
    print!("   -o{}", txt_option_o);
    print!("   -c{}", mytxt_option_c);
    print!("   -b{}", mytxt_option_b);
    print!("   -M{}", txt_option_M);
    print!("   -T{}", txt_option_T);
    print!("   -n{}", mytxt_option_n_est);
    print!("   -aL{}", txt_option_aL);
    print!("   -AL{}", txt_option_AL);
    print!("   -aS{}", txt_option_aS);
    print!("   -AS{}", txt_option_AS);
    print!("   -B{}", txt_option_B);
    print!("   -g{}", txt_option_g);
    print!("   -D{}", mytxt_option_D);
    print!("   -match{}", txt_option_match2);
    print!("   -mismatch{}", txt_option_mismatch2);
    print!("   -gap{}", txt_option_gap2);
    print!("   -gap-ext{}", txt_option_gap_ext);
    print!("   -bak{}", txt_option_bak);
    println!("   -h\tprint this help\n");
    println!("   Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n");
    println!("   If you find cd-hit useful, please kindly cite:\n");
    println!("   {}", cd_hit_ref1);
    println!("   {}", cd_hit_ref2);
    println!("   {}\n\n", cd_hit_ref3);
    std::process::exit(1);
} // END print_usage_454
