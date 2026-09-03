use std::time::{SystemTime, UNIX_EPOCH};

/// Emit the C `__DATE__` / `__TIME__` equivalents that the CD-HIT banner prints.
///
/// Upstream stamps the compiler's build timestamp into the binary.
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (y, m, d, hh, mm, ss) = civil_from_unix(secs);
    const MONTHS: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    // `__DATE__` pads the day with a space, not a zero. `Options::SetOptions` rewrites
    // that space to '0' for its own banner; `cd_hit_ver` keeps the raw form.
    println!(
        "cargo:rustc-env=CDHIT_BUILD_DATE={} {:>2} {}",
        MONTHS[(m - 1) as usize],
        d,
        y
    );
    println!("cargo:rustc-env=CDHIT_BUILD_TIME={:02}:{:02}:{:02}", hh, mm, ss);
}

/// Howard Hinnant's `civil_from_days`, plus the time of day. UTC.
fn civil_from_unix(secs: i64) -> (i64, i64, i64, i64, i64, i64) {
    let days = secs.div_euclid(86400);
    let rem = secs.rem_euclid(86400);
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z.rem_euclid(146097);
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d, rem / 3600, (rem % 3600) / 60, rem % 60)
}
