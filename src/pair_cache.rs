//! Cross-run memory of pairs that cannot pass the diagonal filter at any later threshold.
//!
//! No upstream counterpart. Callers that cluster the *same* sequences repeatedly at falling
//! identity thresholds (Panaroo's `iterative_cdhit`: 0.99, 0.95, ... 0.5, each run on the
//! previous run's representatives) recompute every pair's 2-mer diagonal counts in every
//! run, although those counts do not depend on the threshold. At `-n 2 -c <= 0.8` that
//! diagonal test is ~95% of the time and rejects 95-99.9% of pairs.
//!
//! After a pair's full test, [`RoundCache::mark_if_hopeless`] takes the largest match count
//! over every window any later run could examine. `best_sum` is always the (possibly trimmed)
//! count of one such window, so if that maximum is below the smallest `required_aas` any
//! later run can have for this query, the pair fails every later diagonal test, and one bit
//! records it. Later runs skip those pairs before the test. Pairs are never skipped on
//! anything but that proof, so clustering output is unchanged.
//!
//! Sequences are matched across runs by header identifier. Floors are computed with
//! cd-hit's own `cal_aax_cutoff` + `compute_required_bases`, per query length, assuming
//! `len_eff == len` — which holds because representatives are never shorter than the query;
//! callers check it per query and neither mark nor skip otherwise.

use crate::cdhit_common::{cal_aax_cutoff, Options, WorkingParam, NAA_STAT, NAA_STAT_START_PERCENT};
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

/// Persistent across runs: sequence ids and the "fails every later run" bit per pair.
pub struct PairCache {
    thresholds: Vec<f64>,
    /// the bitset may not grow past this; runs over more sequences go uncached
    max_bytes: usize,
    ids: HashMap<Vec<u8>, u32>,
    /// Lower-triangular bitset over id pairs `(lo < hi)`, bit `hi*(hi-1)/2 + lo`.
    bits: Vec<AtomicU64>,
    band_width: Option<i32>,
    naa: Option<i32>,
    /// Pairs marked (diagnostic). Lookups are deliberately not counted: a shared counter
    /// written on every pair by every thread costs more than the lookup.
    pub marked: AtomicU64,
}

/// One run's view: the run's sequence ids (by `Sequence::index`) and the length floors.
pub struct RoundCache<'a> {
    cache: &'a PairCache,
    /// id of each input sequence, indexed by `Sequence::index`
    ids: Vec<u32>,
    /// no later runs: nothing worth marking
    last: bool,
    /// per query length: min over this and later runs of `required_aa1` / `required_aas`
    floor_aa1: Vec<i32>,
    floor_aas: Vec<i32>,
}

impl PairCache {
    /// `thresholds` is the full sequence of `-c` values the caller will run, in order. The
    /// bitset takes `n * (n - 1) / 16` bytes for `n` distinct sequences, e.g. 13 MB at 15k and
    /// 625 MB at 100k; past `max_bytes` a run gets no cache and behaves as `do_clustering`.
    pub fn new(thresholds: Vec<f64>, max_bytes: usize) -> Self {
        PairCache {
            thresholds,
            max_bytes,
            ids: HashMap::new(),
            bits: Vec::new(),
            band_width: None,
            naa: None,
            marked: AtomicU64::new(0),
        }
    }

    /// Prepare run `round` (an index into the thresholds) over `identifiers` (by
    /// `Sequence::index`), for queries up to `max_len` residues. `None` when the sequences
    /// seen so far would take the bitset past `max_bytes`.
    pub fn begin_round<'a>(
        &'a mut self,
        round: usize,
        identifiers: &[&[u8]],
        max_len: usize,
        options: &Options,
    ) -> Option<RoundCache<'a>> {
        assert!(
            (self.thresholds[round] - options.cluster_thd).abs() < 1e-12,
            "pair cache: run {round} is -c {}, but the cache was set up for {}",
            options.cluster_thd,
            self.thresholds[round]
        );
        // bounds are for this window width and floors for this word length (`-n` enters
        // `compute_required_bases`); runs that change either cannot share marks
        match self.band_width {
            None => self.band_width = Some(options.band_width),
            Some(b) => assert_eq!(b, options.band_width, "pair cache: -b changed between runs"),
        }
        match self.naa {
            None => self.naa = Some(options.naa),
            Some(n) => assert_eq!(n, options.naa, "pair cache: -n changed between runs"),
        }
        let mut ids = Vec::with_capacity(identifiers.len());
        for &name in identifiers {
            let next = self.ids.len() as u32;
            ids.push(*self.ids.entry(name.to_vec()).or_insert(next));
        }
        let n = self.ids.len() as u64;
        let words = (n * n.saturating_sub(1) / 2).div_ceil(64) as usize;
        if words * 8 > self.max_bytes {
            return None;
        }
        if self.bits.len() < words {
            self.bits.resize_with(words, || AtomicU64::new(0));
        }

        // Floors include this run: requirements are not monotone in the threshold (e.g. a
        // 3-residue query needs 0 at c=0.92 but 1 at 0.91), so a mark must also hold for the
        // run that sets it, in case a pair is ever seen twice in one run.
        let from_now = &self.thresholds[round..];
        let last = from_now.len() == 1;
        let (mut floor_aa1, mut floor_aas) = (vec![i32::MAX; max_len + 1], vec![i32::MAX; max_len + 1]);
        for &c in from_now {
            let mut opts = options.clone();
            opts.cluster_thd = c;
            let (mut a1, mut a2, mut an) = (0.0, 0.0, 0.0);
            cal_aax_cutoff(&mut a1, &mut a2, &mut an, c, opts.tolerance, NAA_STAT_START_PERCENT, &NAA_STAT, opts.naa);
            for len in 0..=max_len {
                let mut p = WorkingParam::new(a1, a2, an);
                p.control_short_coverage(len as i32, &opts);
                p.compute_required_bases(opts.naa, 2, &opts);
                floor_aa1[len] = floor_aa1[len].min(p.required_aa1);
                floor_aas[len] = floor_aas[len].min(p.required_aas);
            }
        }
        Some(RoundCache { cache: self, ids, last, floor_aa1, floor_aas })
    }

    fn bit(&self, a: u32, b: u32) -> (usize, u64) {
        let (lo, hi) = if a < b { (a as u64, b as u64) } else { (b as u64, a as u64) };
        let k = hi * (hi - 1) / 2 + lo;
        ((k / 64) as usize, 1u64 << (k % 64))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn required(len: i32, c: f64) -> (i32, i32) {
        let mut opts = Options::default();
        opts.cluster_thd = c;
        opts.naa = 2;
        let (mut a1, mut a2, mut an) = (0.0, 0.0, 0.0);
        cal_aax_cutoff(&mut a1, &mut a2, &mut an, c, opts.tolerance, NAA_STAT_START_PERCENT, &NAA_STAT, 2);
        let mut p = WorkingParam::new(a1, a2, an);
        p.control_short_coverage(len, &opts);
        p.compute_required_bases(2, 2, &opts);
        (p.required_aa1, p.required_aas)
    }

    /// Requirements are not monotone in the threshold, which is why floors include the
    /// marking run itself (see `begin_round`). If this starts failing, that comment is stale.
    #[test]
    fn requirements_can_rise_as_threshold_falls() {
        assert!(required(3, 0.91).1 > required(3, 0.92).1);
    }

    /// Every floor must be at most the run's own requirement, for every length.
    #[test]
    fn floors_never_exceed_any_covered_run() {
        let thresholds = vec![0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5];
        let mut cache = PairCache::new(thresholds.clone(), usize::MAX);
        for round in 0..thresholds.len() {
            let mut opts = Options::default();
            opts.cluster_thd = thresholds[round];
            opts.naa = 2;
            let r = cache.begin_round(round, &[], 3000, &opts).unwrap();
            for len in 0..3000 {
                for &c in &thresholds[round..] {
                    let (aa1, aas) = required(len as i32, c);
                    assert!(r.floor_aa1[len] <= aa1 && r.floor_aas[len] <= aas, "len {len} c {c}");
                }
            }
        }
    }

    #[test]
    fn marks_are_symmetric_and_distinct() {
        let mut cache = PairCache::new(vec![0.9, 0.5], usize::MAX);
        let mut opts = Options::default();
        opts.cluster_thd = 0.9;
        opts.naa = 2;
        let names: Vec<&[u8]> = vec![b">a", b">b", b">c"];
        let r = cache.begin_round(0, &names, 50, &opts).unwrap();
        r.mark_if_hopeless(0, 2, 40, -1);
        assert!(r.is_hopeless(2, 0) && r.is_hopeless(0, 2));
        assert!(!r.is_hopeless(0, 1) && !r.is_hopeless(1, 2));
    }

    #[test]
    fn over_budget_runs_go_uncached() {
        let mut cache = PairCache::new(vec![0.9, 0.5], 0);
        let mut opts = Options::default();
        opts.cluster_thd = 0.9;
        let names: Vec<&[u8]> = vec![b">a", b">b"];
        assert!(cache.begin_round(0, &names, 50, &opts).is_none());
    }
}

impl RoundCache<'_> {
    /// Whether an earlier run proved this pair fails the diagonal test at this threshold.
    pub fn is_hopeless(&self, query_index: i32, rep_index: i32) -> bool {
        let (a, b) = (self.ids[query_index as usize], self.ids[rep_index as usize]);
        if a == b {
            return false;
        }
        let (w, m) = self.cache.bit(a, b);
        self.cache.bits[w].load(Relaxed) & m != 0
    }

    /// Diagonal range `[lo, hi]` this and every later run's band lies inside, for a query of
    /// `len1`.
    ///
    /// One diagonal wider at the low end than the band itself: the bit is shared by both
    /// orientations of a pair, and if a later run tests two equal-length sequences the other
    /// way round, diagonal `d` becomes `nall - 1 - d`, so its band `[b, nall - b]` is
    /// `[b - 1, nall - b - 1]` in this run's numbering.
    pub fn later_band(&self, len1: i32, nall: i32) -> Option<(i32, i32)> {
        if self.last {
            return None;
        }
        let aa1 = *self.floor_aa1.get(len1 as usize)?;
        let lo = (aa1 - 2).max(0);
        (lo <= nall - lo).then_some((lo, nall - lo))
    }

    /// Record the pair as hopeless if `max_window_count` (the largest match count of any
    /// `band_width`-wide window inside [`Self::later_band`]) is below every later floor.
    pub fn mark_if_hopeless(&self, query_index: i32, rep_index: i32, len1: i32, max_window_count: i32) {
        let Some(&floor) = self.floor_aas.get(len1 as usize) else { return };
        if self.last || max_window_count >= floor {
            return;
        }
        let (a, b) = (self.ids[query_index as usize], self.ids[rep_index as usize]);
        if a == b {
            return;
        }
        let (w, m) = self.cache.bit(a, b);
        if self.cache.bits[w].fetch_or(m, Relaxed) & m == 0 {
            self.cache.marked.fetch_add(1, Relaxed);
        }
    }
}
