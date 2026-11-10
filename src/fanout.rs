//! A few persistent threads that run one borrowed closure at a time.
//!
//! No upstream counterpart. `std::thread::scope` costs a thread spawn per worker per call,
//! too much for the few dozen sequences [`crate::window`] hands out at a time; these
//! threads stay up between calls and spin briefly before sleeping, so back-to-back calls
//! start within microseconds.

use std::panic::{catch_unwind, resume_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::*};
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;

/// `(worker, piece)`: `worker` is 0 for the caller and `1..` for the helpers.
type Task = dyn Fn(usize, usize) + Sync;

/// How long a helper spins for the next call before it sleeps.
const SPINS: u32 = 20_000;

struct State {
    task: Option<&'static Task>,
    pieces: usize,
    /// helpers that have not finished the current call
    running: usize,
    panic: Option<Box<dyn std::any::Any + Send>>,
    shutdown: bool,
}

struct Shared {
    state: Mutex<State>,
    wake: Condvar,
    done: Condvar,
    /// bumped once per call; lets a spinning helper see a new call without the lock
    epoch: AtomicU64,
    next_piece: AtomicUsize,
}

pub struct Fanout {
    shared: Arc<Shared>,
    helpers: Vec<JoinHandle<()>>,
}

impl Fanout {
    pub fn new(helpers: usize) -> Self {
        let shared = Arc::new(Shared {
            state: Mutex::new(State { task: None, pieces: 0, running: 0, panic: None, shutdown: false }),
            wake: Condvar::new(),
            done: Condvar::new(),
            epoch: AtomicU64::new(0),
            next_piece: AtomicUsize::new(0),
        });
        let helpers = (1..=helpers)
            .map(|worker| {
                let shared = Arc::clone(&shared);
                std::thread::spawn(move || helper_loop(&shared, worker))
            })
            .collect();
        Fanout { shared, helpers }
    }

    /// Caller plus helpers.
    pub fn workers(&self) -> usize {
        self.helpers.len() + 1
    }

    /// Call `task(worker, piece)` once for each piece in `0..pieces`, spread over the caller
    /// and the helpers; returns once every call has returned. A panic in `task` is re-raised
    /// here.
    pub fn run(&mut self, pieces: usize, task: &(dyn Fn(usize, usize) + Sync)) {
        // SAFETY: the reference is dropped from `State` before this function returns, and it
        // does not return (or unwind) until every helper has finished with it. `&mut self`
        // rules out a second call, from another thread or from inside `task`, resetting
        // `running` while helpers still hold the reference.
        let task: &'static Task = unsafe { std::mem::transmute(task) };
        let shared = &*self.shared;
        {
            let mut state = shared.state.lock().unwrap();
            shared.next_piece.store(0, Relaxed);
            state.task = Some(task);
            state.pieces = pieces;
            state.running = self.helpers.len();
            shared.epoch.fetch_add(1, Release);
        }
        shared.wake.notify_all();
        let own = catch_unwind(AssertUnwindSafe(|| drain(shared, pieces, 0, task)));
        let mut state = shared.state.lock().unwrap();
        while state.running != 0 {
            state = shared.done.wait(state).unwrap();
        }
        state.task = None;
        let panic = state.panic.take();
        drop(state);
        if let Some(payload) = own.err().or(panic) {
            resume_unwind(payload);
        }
    }
}

impl Drop for Fanout {
    fn drop(&mut self) {
        self.shared.state.lock().unwrap().shutdown = true;
        self.shared.epoch.fetch_add(1, Release);
        self.shared.wake.notify_all();
        for helper in self.helpers.drain(..) {
            let _ = helper.join();
        }
    }
}

fn drain(shared: &Shared, pieces: usize, worker: usize, task: &Task) {
    loop {
        let piece = shared.next_piece.fetch_add(1, Relaxed);
        if piece >= pieces {
            break;
        }
        task(worker, piece);
    }
}

fn helper_loop(shared: &Shared, worker: usize) {
    let mut seen = 0u64;
    loop {
        let mut spins = 0;
        while shared.epoch.load(Acquire) == seen && spins < SPINS {
            std::hint::spin_loop();
            spins += 1;
        }
        let (task, pieces) = {
            let mut state = shared.state.lock().unwrap();
            while shared.epoch.load(Acquire) == seen && !state.shutdown {
                state = shared.wake.wait(state).unwrap();
            }
            if state.shutdown {
                return;
            }
            seen = shared.epoch.load(Acquire);
            (state.task.expect("a call is in progress"), state.pieces)
        };
        let outcome = catch_unwind(AssertUnwindSafe(|| drain(shared, pieces, worker, task)));
        let mut state = shared.state.lock().unwrap();
        if let Err(payload) = outcome {
            state.panic.get_or_insert(payload);
        }
        state.running -= 1;
        if state.running == 0 {
            shared.done.notify_one();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_piece_runs_once_across_many_calls() {
        let mut fanout = Fanout::new(3);
        let workers = fanout.workers();
        for pieces in [0usize, 1, 2, 7, 64] {
            for _ in 0..200 {
                let hits: Vec<AtomicUsize> = (0..pieces).map(|_| AtomicUsize::new(0)).collect();
                fanout.run(pieces, &|worker, piece| {
                    assert!(worker < workers);
                    hits[piece].fetch_add(1, Relaxed);
                });
                assert!(hits.iter().all(|h| h.load(Relaxed) == 1));
            }
        }
    }

    #[test]
    fn a_panicking_piece_is_reported_and_the_pool_survives() {
        let mut fanout = Fanout::new(2);
        let caught = catch_unwind(AssertUnwindSafe(|| {
            fanout.run(8, &|_, piece| assert_ne!(piece, 5, "boom"));
        }));
        assert!(caught.is_err());
        let count = AtomicUsize::new(0);
        fanout.run(8, &|_, _| {
            count.fetch_add(1, Relaxed);
        });
        assert_eq!(count.load(Relaxed), 8);
    }
}
