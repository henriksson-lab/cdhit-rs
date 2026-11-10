//! The representative's packed 4-mer codes for the nucleotide diagonal kernel.
//!
//! Position `j` (of `iseq2.len() - 3`) gets `c0<<6 | c1<<4 | c2<<2 | c3` from `iseq2[j..j+4]`
//! when all four are bases (0..=3); otherwise `ra` keeps its 0xFF and `rb` its 0x00 sentinel.
//! The vector form writes the sentinels itself, which is the same since the caller keeps both
//! buffers at their sentinels outside the encoded span.

/// Encode into `ra`/`rb` from offset `at`. `simd` must only be true with AVX-512F and BW.
pub(crate) fn encode(iseq2: &[u8], ra: &mut [u8], rb: &mut [u8], at: usize, simd: bool) {
    let nr = iseq2.len().saturating_sub(3);
    assert!(at + nr <= ra.len() && at + nr <= rb.len());
    let mut j = 0;
    #[cfg(target_arch = "x86_64")]
    if simd {
        // SAFETY: AVX-512F/BW per the caller; loads end at `j + 3 + 63 < nr + 3 = iseq2.len()`
        // and stores at `at + j + 63 < at + nr`, both checked above.
        j = unsafe { encode_avx512(iseq2, ra, rb, at, nr) };
    }
    let _ = simd;
    encode_scalar(iseq2, ra, rb, at, j, nr);
}

fn encode_scalar(iseq2: &[u8], ra: &mut [u8], rb: &mut [u8], at: usize, from: usize, nr: usize) {
    for (j, w) in iseq2.windows(4).enumerate().take(nr).skip(from) {
        if (w[0] | w[1] | w[2] | w[3]) >= 4 {
            continue;
        }
        let code = (w[0] << 6) | (w[1] << 4) | (w[2] << 2) | w[3];
        ra[at + j] = code;
        rb[at + j] = code;
    }
}

/// Whole 64-position blocks; returns the first position not encoded.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn encode_avx512(iseq2: &[u8], ra: &mut [u8], rb: &mut [u8], at: usize, nr: usize) -> usize {
    use core::arch::x86_64::*;
    let p = iseq2.as_ptr();
    let three = _mm512_set1_epi8(3);
    let ff = _mm512_set1_epi8(-1);
    let mut j = 0;
    while j + 64 <= nr {
        let c: [__m512i; 4] = core::array::from_fn(|k| _mm512_loadu_si512(p.add(j + k) as *const _));
        let bad = _mm512_or_si512(_mm512_or_si512(c[0], c[1]), _mm512_or_si512(c[2], c[3]));
        let valid = _mm512_testn_epi8_mask(bad, _mm512_set1_epi8(!3));
        // Masked to 2 bits first so the 16-bit shifts cannot carry into the neighbouring byte.
        let s = |k: usize, by: u32| match by {
            6 => _mm512_slli_epi16::<6>(_mm512_and_si512(c[k], three)),
            4 => _mm512_slli_epi16::<4>(_mm512_and_si512(c[k], three)),
            2 => _mm512_slli_epi16::<2>(_mm512_and_si512(c[k], three)),
            _ => _mm512_and_si512(c[k], three),
        };
        let code = _mm512_or_si512(_mm512_or_si512(s(0, 6), s(1, 4)), _mm512_or_si512(s(2, 2), s(3, 0)));
        _mm512_storeu_si512(ra.as_mut_ptr().add(at + j) as *mut _, _mm512_mask_blend_epi8(valid, ff, code));
        _mm512_storeu_si512(rb.as_mut_ptr().add(at + j) as *mut _, _mm512_maskz_mov_epi8(valid, code));
        j += 64;
    }
    j
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simd_matches_scalar() {
        if !(cfg!(target_arch = "x86_64")
            && std::is_x86_feature_detected!("avx512f")
            && std::is_x86_feature_detected!("avx512bw"))
        {
            return;
        }
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = |m: u64| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s % m
        };
        for case in 0..3000 {
            let n = next(400) as usize;
            // mostly bases, some N (4) and mask (5) codes, and the odd high byte
            let seq: Vec<u8> = (0..n)
                .map(|_| match next(100) {
                    0..=89 => next(4) as u8,
                    90..=95 => 4,
                    96..=98 => 5,
                    _ => next(256) as u8,
                })
                .collect();
            let at = next(70) as usize;
            let len = at + n + 8;
            let (mut a0, mut b0) = (vec![0xFFu8; len], vec![0u8; len]);
            let (mut a1, mut b1) = (a0.clone(), b0.clone());
            encode(&seq, &mut a0, &mut b0, at, false);
            encode(&seq, &mut a1, &mut b1, at, true);
            assert_eq!((a0, b0), (a1, b1), "case {case} n={n} at={at}");
        }
    }
}
