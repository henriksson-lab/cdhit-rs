//! The sliding band-window scan at the end of the diagonal tests. Not upstream.
//!
//! Upstream slides a `w`-wide window across the per-diagonal counts one diagonal at a time,
//! recording the window whenever its complexity-weighted sum (`score2`) strictly beats the
//! best so far. That is a serial add/compare chain per diagonal, yet a new best is rare once
//! the window has passed the main hit. On AVX-512 the scan computes 16 consecutive window sums
//! with an in-register prefix sum and resolves the block's records (see `slide_avx512`)
//! without the per-diagonal dependency; ties, `imax` updates and every output match upstream.

/// Running state of the window scan in `finish_diag_band_u16`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Slide {
    pub score: i32,
    pub score2: i32,
    pub best_score: i32,
    pub best_score2: i32,
    pub max_diag2: i32,
    pub imax_off: usize,
    pub from_off: usize,
    pub end_off: usize,
}

/// Slide the `w`-wide window so its right edge visits diagonals `j0..j1`; the left edge
/// drops diagonal `j - w`. `simd` must only be true when AVX-512F is available.
pub(crate) fn slide(cnt: &[u16], cpx: &[u16], w: usize, j0: usize, j1: usize, st: &mut Slide, simd: bool) {
    assert!(w <= j0 && j0 <= j1 && j1 <= cnt.len() && j1 <= cpx.len());
    #[cfg(target_arch = "x86_64")]
    if simd {
        // SAFETY: the caller guarantees AVX-512F; the bounds were asserted above.
        unsafe { slide_avx512(cnt, cpx, w, j0, j1, st) };
        return;
    }
    let _ = simd;
    slide_scalar(cnt, cpx, w, j0, j1, st);
}

fn slide_scalar(cnt: &[u16], cpx: &[u16], w: usize, j0: usize, j1: usize, st: &mut Slide) {
    let entering = cnt[j0..j1].iter().zip(&cpx[j0..j1]);
    let leaving = cnt[j0 - w..j1 - w].iter().zip(&cpx[j0 - w..j1 - w]);
    for (j, ((&jc, &jx), (&kc, &kx))) in (j0..).zip(entering.zip(leaving)) {
        st.score += jc as i32 - kc as i32;
        st.score2 += jx as i32 - kx as i32;
        if st.score2 > st.best_score2 {
            st.from_off = j - w + 1;
            st.end_off = j;
            st.best_score = st.score;
            st.best_score2 = st.score2;
            if jx as i32 > st.max_diag2 {
                st.max_diag2 = jx as i32;
                st.imax_off = j;
            }
        }
    }
}

/// Inclusive prefix sum across the 16 lanes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn prefix_sum(v: core::arch::x86_64::__m512i) -> core::arch::x86_64::__m512i {
    use core::arch::x86_64::*;
    let z = _mm512_setzero_si512();
    let mut p = v;
    p = _mm512_add_epi32(p, _mm512_alignr_epi32::<15>(p, z));
    p = _mm512_add_epi32(p, _mm512_alignr_epi32::<14>(p, z));
    p = _mm512_add_epi32(p, _mm512_alignr_epi32::<12>(p, z));
    _mm512_add_epi32(p, _mm512_alignr_epi32::<8>(p, z))
}

/// One block = 16 scalar steps at once. Lane `l` is a record when its window sum beats both
/// `best_score2` and every earlier lane; the last record fixes the band, and the scalar
/// `if cpx > max_diag2` chain over the records ends on the first record lane holding their
/// largest `cpx`, so `imax` is that lane if its value beats `max_diag2`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn slide_avx512(cnt: &[u16], cpx: &[u16], w: usize, j0: usize, j1: usize, st: &mut Slide) {
    use core::arch::x86_64::*;
    let widen = |p: &[u16], at: usize| _mm512_cvtepu16_epi32(_mm256_loadu_si256(p.as_ptr().add(at) as *const __m256i));
    let floor = _mm512_set1_epi32(i32::MIN);
    let mut j = j0;
    while j + 16 <= j1 {
        let k = j - w;
        let x2 = widen(cpx, j);
        let d2 = _mm512_sub_epi32(x2, widen(cpx, k));
        let s2 = _mm512_add_epi32(prefix_sum(d2), _mm512_set1_epi32(st.score2));
        let best2 = _mm512_set1_epi32(st.best_score2);
        let mut pm = s2;
        pm = _mm512_max_epi32(pm, _mm512_alignr_epi32::<15>(pm, floor));
        pm = _mm512_max_epi32(pm, _mm512_alignr_epi32::<14>(pm, floor));
        pm = _mm512_max_epi32(pm, _mm512_alignr_epi32::<12>(pm, floor));
        pm = _mm512_max_epi32(pm, _mm512_alignr_epi32::<8>(pm, floor));
        let before = _mm512_max_epi32(best2, _mm512_alignr_epi32::<15>(pm, best2));
        let rec = _mm512_cmpgt_epi32_mask(s2, before);
        let d = _mm512_sub_epi32(widen(cnt, j), widen(cnt, k));
        if rec == 0 {
            st.score2 += _mm512_reduce_add_epi32(d2);
            st.score += _mm512_reduce_add_epi32(d);
        } else {
            let mut s2v = [0i32; 16];
            let mut sv = [0i32; 16];
            _mm512_storeu_si512(s2v.as_mut_ptr() as *mut __m512i, s2);
            _mm512_storeu_si512(sv.as_mut_ptr() as *mut __m512i, _mm512_add_epi32(prefix_sum(d), _mm512_set1_epi32(st.score)));
            let last = 15 - rec.leading_zeros() as usize;
            st.best_score2 = s2v[last];
            st.best_score = sv[last];
            st.end_off = j + last;
            st.from_off = j + last - w + 1;
            let mx = _mm512_mask_reduce_max_epi32(rec, x2);
            if mx > st.max_diag2 {
                let at = _mm512_mask_cmpeq_epi32_mask(rec, x2, _mm512_set1_epi32(mx));
                st.max_diag2 = mx;
                st.imax_off = j + at.trailing_zeros() as usize;
            }
            st.score2 = s2v[15];
            st.score = sv[15];
        }
        j += 16;
    }
    slide_scalar(cnt, cpx, w, j, j1, st);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn xorshift(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }

    fn start(cnt: &[u16], cpx: &[u16], w: usize) -> Slide {
        let score = cnt[..w].iter().map(|&c| c as i32).sum();
        let score2 = cpx[..w].iter().map(|&c| c as i32).sum();
        let (imax_off, max_diag2) = cpx[..w]
            .iter()
            .enumerate()
            .fold((0, 0), |(bi, bv), (i, &v)| if (v as i32) > bv { (i, v as i32) } else { (bi, bv) });
        Slide { score, score2, best_score: score, best_score2: score2, max_diag2, imax_off, from_off: 0, end_off: w - 1 }
    }

    /// Sparse counts with occasional peaks and plateaus, like real diagonals, so blocks both
    /// skip and fall back; plateaus exercise the strict `>` tie rule.
    #[test]
    fn simd_matches_scalar() {
        if !(cfg!(target_arch = "x86_64") && std::is_x86_feature_detected!("avx512f")) {
            return;
        }
        let mut s = 0x9E37_79B9_7F4A_7C15u64;
        for case in 0..20_000 {
            let n = 2 + (xorshift(&mut s) % 900) as usize;
            let mut cnt = vec![0u16; n];
            let mut cpx = vec![0u16; n];
            let peak = (xorshift(&mut s) % n as u64) as usize;
            for i in 0..n {
                let r = xorshift(&mut s);
                let near = peak.abs_diff(i) < 1 + (r % 30) as usize;
                let c = match r % 7 {
                    _ if near => (r >> 8) % 200,
                    0 | 1 => (r >> 8) % 3,
                    2 if case % 3 == 0 => 5,
                    _ => 0,
                } as u16;
                cnt[i] = c;
                cpx[i] = c + if case % 2 == 0 { c } else { ((r >> 20) % (c as u64 + 1)) as u16 };
            }
            let w = 1 + (xorshift(&mut s) % (n as u64 - 1).min(40)) as usize;
            let j1 = n - (xorshift(&mut s) % 2) as usize;
            if j1 < w {
                continue;
            }
            let mut a = start(&cnt, &cpx, w);
            let mut b = a;
            slide(&cnt, &cpx, w, w, j1, &mut a, false);
            slide(&cnt, &cpx, w, w, j1, &mut b, true);
            assert_eq!(a, b, "case {case} n={n} w={w} j1={j1}");
        }
    }
}
