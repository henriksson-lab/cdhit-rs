//! Dense word counts for protein `-n 2` tables. Not upstream.
//!
//! At `-n 2` there are only 441 words and nearly every representative has most of them, so
//! `count_words_direct`'s posting-list walk visits about one entry per (representative, query
//! word) -- the family pass's single most expensive loop. The shared-word count it computes is
//! `sum_w min(query[w], rep[w])`, which with one byte per word is seven AVX-512 `min` + `sad`
//! steps per representative.
//!
//! The walk's result is reproduced exactly:
//! * Candidates are the representatives whose count reaches `min`. The walk only admits a
//!   representative while `rest = aan_no - k + 1 >= min` at its first shared word `k`, but one
//!   it refuses can gain at most `aan_no - k < min` from the words left, so it would fail the
//!   final `count >= min` filter anyway.
//! * Candidates are ordered by that first shared word, then by table index. The query's words
//!   are sorted, so the first shared word is the representative's lowest shared code.
//! * Bytes saturate at 255. A representative's saturated count still gives the exact `min`
//!   against any query count up to 255; queries with a larger count use the walk.

/// Words at `-n 2` over cd-hit's 21-letter protein alphabet.
pub(crate) const CODES: usize = 441;
const ROW: usize = 448;

/// One byte row per representative, by table index.
#[derive(Default)]
pub(crate) struct DenseWords {
    rows: Vec<u8>,
    n: usize,
    live: bool,
}

/// What [`DenseWords::count`] needs of the query: its sorted word codes and their run
/// multiplicities, as `encode_words` leaves them.
pub(crate) struct Query<'a> {
    pub words: &'a [i32],
    pub counts: &'a [u32],
}

impl DenseWords {
    /// Empty, and live only for a table of [`CODES`] words.
    pub fn reset(&mut self, naan: i32) {
        self.rows.clear();
        self.n = 0;
        self.live = naan as usize == CODES;
    }

    pub fn is_live(&self) -> bool {
        self.live
    }

    /// Add the representative at table index `idx` with these `(word, count)` entries. Any
    /// layout the walk would count differently (an out-of-order index, a word outside the
    /// alphabet, or a word listed twice) turns the dense counts off for this table.
    pub fn add(&mut self, idx: i32, entries: impl Iterator<Item = (i32, u32)>) {
        if !self.live {
            return;
        }
        if idx as usize != self.n {
            return self.kill();
        }
        let base = self.rows.len();
        self.rows.resize(base + ROW, 0);
        for (word, count) in entries {
            if word < 0 || word as usize >= CODES || self.rows[base + word as usize] != 0 || count == 0 {
                return self.kill();
            }
            self.rows[base + word as usize] = count.min(255) as u8;
        }
        self.n += 1;
    }

    fn kill(&mut self) {
        self.live = false;
        self.rows = Vec::new();
    }

    /// The candidates `count_words_direct` would list for this query among table indices
    /// `from..`, as `(index, count)` in its order; `None` if the query needs the walk.
    /// `simd` must only be true when AVX-512F and BW are available.
    pub fn count(&self, q: &Query, from: usize, min: i32, simd: bool) -> Option<Vec<(i32, i32)>> {
        if !self.live || !simd {
            return None;
        }
        let mut cq = [0u8; ROW];
        let mut first_k = [0u32; CODES];
        let mut prev = -1i32;
        for (k, (&w, &c)) in q.words.iter().zip(q.counts).enumerate() {
            if c == 0 {
                continue;
            }
            if w <= prev || w as usize >= CODES || c > 255 {
                return None;
            }
            prev = w;
            cq[w as usize] = c as u8;
            first_k[w as usize] = k as u32;
        }
        let threshold = min.max(1) as u64;
        #[cfg(target_arch = "x86_64")]
        {
            let mut hits: Vec<(u64, i32)> = Vec::new();
            // SAFETY: AVX-512F/BW per the caller; rows `from..n` are inside `self.rows`.
            unsafe { scan_avx512(&self.rows, from, self.n, &cq, threshold, &mut hits) };
            let mut out: Vec<(u64, i32)> = hits
                .into_iter()
                .map(|(hit, s)| ((first_k[(hit >> 32) as usize] as u64) << 32 | (hit & 0xffff_ffff), s))
                .collect();
            out.sort_unstable_by_key(|&(key, _)| key);
            Some(out.into_iter().map(|(key, s)| (key as u32 as i32, s)).collect())
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (cq, first_k, threshold, from);
            None
        }
    }
}

/// Push `(lowest shared code << 32 | index, count)` for every row in `from..n` whose
/// `sum min(row, cq)` reaches `threshold`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn scan_avx512(rows: &[u8], from: usize, n: usize, cq: &[u8; ROW], threshold: u64, out: &mut Vec<(u64, i32)>) {
    use core::arch::x86_64::*;
    debug_assert!(n * ROW <= rows.len());
    let zero = _mm512_setzero_si512();
    let q: [__m512i; 7] = std::array::from_fn(|v| _mm512_loadu_si512(cq.as_ptr().add(v * 64) as *const __m512i));
    let qmask: [u64; 7] = std::array::from_fn(|v| _mm512_test_epi8_mask(q[v], q[v]));
    for r in from..n {
        let p = rows.as_ptr().add(r * ROW);
        let mut acc = zero;
        for (v, &qv) in q.iter().enumerate() {
            let x = _mm512_loadu_si512(p.add(v * 64) as *const __m512i);
            acc = _mm512_add_epi64(acc, _mm512_sad_epu8(_mm512_min_epu8(x, qv), zero));
        }
        let s = _mm512_reduce_add_epi64(acc) as u64;
        if s < threshold {
            continue;
        }
        let code = (0..7)
            .find_map(|v| {
                let x = _mm512_loadu_si512(p.add(v * 64) as *const __m512i);
                let m = _mm512_test_epi8_mask(x, x) & qmask[v];
                (m != 0).then(|| v * 64 + m.trailing_zeros() as usize)
            })
            .expect("a positive shared count has a shared word");
        out.push(((code as u64) << 32 | r as u64, s as i32));
    }
}
