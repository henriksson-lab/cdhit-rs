//! Check a window of a block's sequences at once, then settle them in order.
//!
//! No upstream counterpart. Building a block's word table is one sequence at a time
//! upstream: each sequence must see the representatives made from the ones before it. At
//! `-T 10` on pangenome input that serial phase is close to half the run.
//!
//! # Why the result is unchanged
//!
//! Take the next [`WINDOW`] sequences. The table as it stands holds every representative
//! before the window, so checking a window sequence against it now gives what its own turn
//! would give, *unless* a representative created earlier in the same window would have been
//! accepted for it. Rejected candidates never change a sequence's state, so that is the only
//! way the turn-by-turn result can differ. Hence, per window:
//!
//! 1. every sequence is checked against the table as it stands (all threads);
//! 2. every sequence is checked against the window's would-be representatives that precede
//!    it, starting from a blank state, which accepts exactly the ones that pass the
//!    thresholds on their own (all threads);
//! 3. in order, on one thread: a sequence with no such collision takes its step-1 result,
//!    becoming a representative or not just as `cluster_one` would decide. One with a
//!    collision is put back as it was and goes through `cluster_one` itself, as does one
//!    whose identical twin (see `twin`) is in the window.
//!
//! Step 2's set is every sequence step 1 left without a hit. That covers every
//! representative step 3 can create except one from the `cluster_one` path that step 1 did
//! not predict; after such a one, the rest of the window takes that path too.
//!
//! Step 2 judges each pair as the real table would: [`TwinSkip::before`] hides the window
//! members after the sequence, including from the table-wide lengths that set its
//! thresholds.

use crate::cdhit_common::*;
use crate::fanout::Fanout;
use crate::twin::TwinSkip;
use std::sync::Mutex;

/// Sequences checked per round of helpers. Larger means fewer rounds but more collisions.
const WINDOW: usize = 64;

/// Whether a run's options are ones this reproduces: whole sequences compared by identity,
/// proteins under `-g 1` (best hit) and nucleotides under `-g 0` (first hit).
pub fn covers(options: &Options) -> bool {
    let mode = if options.is_est { !options.cluster_best && !options.is_454 } else { options.cluster_best };
    mode && options.frag_size == 0 && !options.has_2d && !options.use_distance
}

/// Whether to cluster blocks in windows at all. `-aL` is left out because `check_one` can
/// then return before encoding a sequence's words, which step 3 needs.
pub fn worthwhile(options: &Options, threads: usize) -> bool {
    threads > 1 && covers(options) && options.long_coverage == 0.0
}

/// What a block's windows work with.
pub(crate) struct Block<'a, 'g> {
    pub table: &'a mut WordTable,
    /// empty between windows; holds one window's would-be representatives
    pub window_table: &'a mut WordTable,
    /// one per thread
    pub params: &'a mut [WorkingParam],
    pub buffers: &'a mut [WorkingBuffer],
    pub g: Globals<'g>,
    pub twins: &'a [u32],
    pub start: usize,
    /// the table is full at this many entries, or this many representatives
    pub max_items: usize,
    pub max_seqs: usize,
    pub helpers: &'a mut Fanout,
}

/// The parts of a sequence `check_one` can change.
struct Saved {
    state: i16,
    identity: f32,
    cluster_id: i32,
    distance: f32,
    coverage: [i32; 4],
    /// under `-g 0` a hit also frees the residues
    data: Option<Vec<u8>>,
}

impl Saved {
    fn of(seq: &Sequence, options: &Options) -> Self {
        Saved {
            state: seq.state,
            identity: seq.identity,
            cluster_id: seq.cluster_id,
            distance: seq.distance,
            coverage: seq.coverage,
            data: (!options.cluster_best).then(|| seq.data.clone()),
        }
    }

    fn restore(self, seq: &mut Sequence) {
        seq.state = self.state;
        seq.identity = self.identity;
        seq.cluster_id = self.cluster_id;
        seq.distance = self.distance;
        seq.coverage = self.coverage;
        if let Some(data) = self.data {
            seq.data = data;
        }
    }
}

enum Plan {
    /// already redundant: the block loop passes over it
    Skip,
    /// its twin is in the window, so it waits for `cluster_one`
    Serial,
    Checked { view: TwinSkip, saved: Saved },
}

struct Slot {
    plan: Plan,
    /// `check_one`'s return value from step 1
    flag: i32,
    /// its sorted words and their counts, if step 1 left it to become a representative
    words: Option<(Vec<i32>, Vec<INTs>)>,
    collides: bool,
}

type Worker<'a> = Mutex<(&'a mut WorkingParam, &'a mut WorkingBuffer)>;

impl SequenceDB {
    /// The block loop of `do_clustering` over `from..to`, in windows. Returns where the loop
    /// stopped: `to`, or just past the sequence that filled the table.
    pub(crate) fn cluster_block_in_windows(&mut self, block: &mut Block, from: usize, to: usize) -> usize {
        let mut start = from;
        while start < to {
            let end = (start + WINDOW).min(to);
            if let Some(stopped) = self.cluster_window(block, start, end) {
                return stopped;
            }
            start = end;
        }
        to
    }

    /// One window. `Some(i)` if the table filled up at sequence `i - 1`.
    fn cluster_window(&mut self, block: &mut Block, start: usize, end: usize) -> Option<usize> {
        let options = block.g.options;
        let mut slots: Vec<Slot> = (start..end)
            .map(|ks| {
                let plan = if self.sequences[ks].state & IS_REDUNDANT != 0 {
                    Plan::Skip
                } else {
                    match self.twin_skip(block.twins, block.start, start, ks) {
                        None => Plan::Serial,
                        Some(view) => Plan::Checked { view, saved: Saved::of(&self.sequences[ks], options) },
                    }
                };
                Slot { plan, flag: 0, words: None, collides: false }
            })
            .collect();

        self.check_against_table(block, start, &mut slots);
        self.find_collisions(block, start, &mut slots);

        let mut unpredicted_rep = false;
        for k in 0..slots.len() {
            let ks = start + k;
            let slot = &mut slots[k];
            let plan = std::mem::replace(&mut slot.plan, Plan::Skip);
            match plan {
                Plan::Skip => continue,
                Plan::Checked { .. } if !unpredicted_rep && !slot.collides => {
                    self.settle(ks, slot.flag, slot.words.take(), block.table, options);
                }
                plan => {
                    if let Plan::Checked { saved, .. } = plan {
                        saved.restore(&mut self.sequences[ks]);
                    }
                    let reps = self.rep_seqs.len();
                    block.params[0].twin = self.twin_skip(block.twins, block.start, ks, ks).unwrap_or_default();
                    self.cluster_one(ks as i32, block.table, &mut block.params[0], &mut block.buffers[0], block.g);
                    block.params[0].twin = TwinSkip::default();
                    unpredicted_rep |= self.rep_seqs.len() > reps && slot.words.is_none();
                }
            }
            if block.table.size >= block.max_items || block.table.sequences.len() >= block.max_seqs {
                // the rest leave the block untouched, to be screened against its table
                for (later, slot) in slots.iter_mut().enumerate().skip(k + 1) {
                    if let Plan::Checked { saved, .. } = std::mem::replace(&mut slot.plan, Plan::Skip) {
                        saved.restore(&mut self.sequences[start + later]);
                    }
                }
                return Some(ks + 1);
            }
        }
        None
    }

    /// Step 1.
    fn check_against_table(&mut self, block: &mut Block, start: usize, slots: &mut [Slot]) {
        let options = block.g.options;
        let (g, table) = (block.g, &*block.table);
        let (head, tail) = self.sequences.split_at_mut(start);
        let head: &[Sequence] = head;
        let cells = SeqCells::new(&mut tail[..slots.len()]);
        let workers: Vec<Worker> = block.params.iter_mut().zip(block.buffers.iter_mut()).map(Mutex::new).collect();
        let slots: Vec<Mutex<&mut Slot>> = slots.iter_mut().map(Mutex::new).collect();
        block.helpers.run(slots.len(), &|worker, k| {
            let mut slot = slots[k].lock().unwrap();
            let Plan::Checked { view, .. } = slot.plan else { return };
            // SAFETY: `Fanout::run` passes each `k` to one call only.
            let seq = unsafe { cells.get(k) };
            let mut worker = workers[worker].lock().unwrap();
            let (param, buffer) = &mut *worker;
            param.twin = view;
            slot.flag = check_one(seq, head, table, param, buffer, g);
            param.twin = TwinSkip::default();
            // `cluster_one`'s condition for a new representative
            if slot.flag == 0 && !(seq.identity > 0.0 && options.cluster_best) {
                let words = (seq.size - options.naa + 1).max(0) as usize;
                slot.words = Some((buffer.word_encodes[..words].to_vec(), buffer.word_encodes_no[..words].to_vec()));
            }
        });
    }

    /// Step 2.
    fn find_collisions(&mut self, block: &mut Block, start: usize, slots: &mut [Slot]) {
        let options = block.g.options;
        let window_table = &mut *block.window_table;
        for (k, slot) in slots.iter().enumerate() {
            if let Some((words, counts)) = &slot.words {
                let index = window_table.sequences.len() as i32;
                window_table.add_word_counts(words.len() as i32, words, counts, index, options.is_est);
                window_table.sequences.push(start + k);
            }
        }
        let Some(&first) = window_table.sequences.first() else { return };
        {
            let (g, window_table) = (block.g, &*window_table);
            let db: &[Sequence] = &self.sequences;
            let workers: Vec<Worker> =
                block.params.iter_mut().zip(block.buffers.iter_mut()).map(Mutex::new).collect();
            let slots: Vec<Mutex<&mut Slot>> = slots.iter_mut().map(Mutex::new).collect();
            block.helpers.run(slots.len(), &|worker, k| {
                let ks = start + k;
                let mut slot = slots[k].lock().unwrap();
                let Plan::Checked { saved, .. } = &slot.plan else { return };
                if ks <= first {
                    return;
                }
                let real = &db[ks];
                let mut probe = Sequence {
                    data: saved.data.clone().unwrap_or_else(|| real.data.clone()),
                    size: real.size,
                    index: real.index,
                    ..Default::default()
                };
                let mut worker = workers[worker].lock().unwrap();
                let (param, buffer) = &mut *worker;
                param.twin = TwinSkip { before: ks, ..Default::default() };
                slot.collides = check_one(&mut probe, db, window_table, param, buffer, g) != 0;
                param.twin = TwinSkip::default();
            });
        }
        let words = slots.iter().filter_map(|slot| slot.words.as_ref()).flat_map(|(words, _)| words.iter().copied());
        window_table.clear_words(words);
    }

    /// What `cluster_one` does with `check_one`'s result.
    fn settle(
        &mut self,
        ks: usize,
        flag: i32,
        words: Option<(Vec<i32>, Vec<INTs>)>,
        table: &mut WordTable,
        options: &Options,
    ) {
        let seq = &mut self.sequences[ks];
        if flag == 0 {
            if seq.identity > 0.0 && options.cluster_best {
                seq.state |= IS_REDUNDANT;
                seq.clear();
            } else {
                let (words, counts) = words.expect("step 1 kept the words of every would-be representative");
                seq.cluster_id = self.rep_seqs.len() as i32;
                seq.identity = 0.0;
                seq.state |= IS_REP;
                self.rep_seqs.push(ks as i32);
                table.add_word_counts(words.len() as i32, &words, &counts, table.sequences.len() as i32, options.is_est);
                table.sequences.push(ks);
            }
        }
        self.report_progress(ks as i32);
    }
}
