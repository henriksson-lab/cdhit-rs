//! Let a sequence identical to an earlier one start from that sequence's result.
//!
//! No upstream counterpart. Under `-g 1` every sequence is tested against every
//! representative in its length window, and the outcome of one test depends only on the two
//! sequences. So for a sequence `j` whose data equals that of an earlier sequence `t` (its
//! *twin*), every test against a representative before `t` repeats a test `t` already ran.
//! On pangenome input most sequences are such copies.
//!
//! # Why the result is unchanged
//!
//! A sequence's final state is a fold over its candidates in evaluation order (tables in
//! block order, then candidate-list order inside a table), keeping the first candidate of
//! strictly greatest identity. Candidate-list order is (first query word that hits the
//! representative, table index), so two identical queries see their common candidates in the
//! same relative order, and `t`'s final state is the fold over the representatives before
//! `t`. `j` therefore:
//!
//! * is not screened at all while `t` is still unfinished (same tables, same answers);
//! * takes `t`'s state once `t` is final, and from then on skips candidates before `t`.
//!
//! Skipping is done in the word count itself (`count_words_direct`'s `from`), so the skipped
//! representatives cost nothing.
//!
//! That leaves ties. Folding from `t`'s best `O` at identity `x` keeps `O` when a later
//! representative also scores `x`, which is wrong exactly when the real order meets that
//! representative first: `O` sits in the table being read ([`TwinSkip::best_in_table`]) and
//! shares its first word with the query later than the newcomer does. `check_one_aa`
//! compares the two first-word positions when such a tie occurs. Every candidate before `O`
//! scores below `x` (else it would be `t`'s best), so no other tie can go differently.
//!
//! The same reasoning gives the return flag of `check_one` ("accepted a candidate from this
//! table"): the real run accepts one iff a skipped-over `O` is in the table or a newcomer
//! wins.

use crate::cdhit_common::{Options, Sequence, WordTable};
use std::collections::HashMap;

pub const NO_TWIN: u32 = u32::MAX;
const NO_REP: usize = usize::MAX;

/// Which of a table's representatives one `check_one` call looks at; the default is all.
#[derive(Clone, Copy)]
pub struct TwinSkip {
    /// candidates at database positions below this were judged for the twin
    pub below: usize,
    /// database position of the representative the inherited state points at
    pub best_rep: usize,
    /// Candidates at database positions from this on are treated as absent; see `window`.
    pub before: usize,
}

impl Default for TwinSkip {
    fn default() -> Self {
        TwinSkip { below: 0, best_rep: NO_REP, before: usize::MAX }
    }
}

impl TwinSkip {
    /// How many of `table`'s representatives, which are in position order, lie before
    /// [`Self::before`].
    pub fn table_len(&self, table: &WordTable) -> usize {
        if self.before == usize::MAX {
            return table.sequences.len();
        }
        table.sequences.partition_point(|&r| r < self.before)
    }

    /// Whether the inherited best is one of `table_first..`, the positions held by the table
    /// being read.
    pub fn best_in_table(&self, table_first: Option<usize>) -> bool {
        self.best_rep != NO_REP && table_first.is_some_and(|first| self.best_rep >= first)
    }
}

/// Whether `-g 1` is on in the form where a hit at identity 1 can never be displaced:
/// global identity is matches over the query's length, so it cannot exceed 1.
pub fn best_hit_is_capped(options: &Options) -> bool {
    options.cluster_best && options.global_identity && !options.use_distance && !options.has_2d
}

/// For each sequence, the nearest earlier sequence with identical data, or [`NO_TWIN`].
/// Empty when the run's options fall outside what the module argument covers.
pub fn find_twins(sequences: &[Sequence], options: &Options) -> Vec<u32> {
    let covered = best_hit_is_capped(options)
        && !options.is_est
        && options.frag_size == 0
        // `check_one_aa` drops a whole table on `-aL` using the table's shortest
        // representative, which differs between a sequence and its twin
        && options.long_coverage == 0.0
        && sequences.len() < NO_TWIN as usize
        && sequences.windows(2).all(|w| w[0].size >= w[1].size);
    if !covered {
        return Vec::new();
    }
    let mut last: HashMap<&[u8], u32> = HashMap::with_capacity(sequences.len());
    sequences
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            if seq.state != 0 || seq.data.is_empty() {
                return NO_TWIN;
            }
            last.insert(&seq.data, i as u32).unwrap_or(NO_TWIN)
        })
        .collect()
}

/// Copy the twin's result into `seq` and say what `seq` may now skip.
pub fn inherit(seq: &mut Sequence, twin_pos: usize, twin: &Sequence, rep_seqs: &[i32]) -> TwinSkip {
    seq.identity = twin.identity;
    seq.cluster_id = twin.cluster_id;
    seq.distance = twin.distance;
    seq.coverage = twin.coverage;
    let best_rep = if twin.identity > 0.0 { rep_seqs[twin.cluster_id as usize] as usize } else { NO_REP };
    TwinSkip { below: twin_pos, best_rep, ..Default::default() }
}
