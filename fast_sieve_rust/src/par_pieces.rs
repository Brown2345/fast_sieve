#![allow(unused_imports)]
use crate::dioph_appr::{dioph_appr, Rat};
use crate::simple_sieve::simple_seg_sieve;
use core::f64;
use rayon::prelude::*;
use rug::{Integer, Rational};

pub fn s_piece(miter: i64, n: usize, delta: usize) -> Vec<usize> {
    //sieves on a small interval to parallelize code
    let mut _s = vec![1; 2 * delta + 1];
    if miter % 2 == 1 {
        let np = ((n + delta) as f64 / miter as f64).floor() as i64 * miter;
        if np >= (n - delta) as i64 && np <= (n + delta) as i64 && np > miter {
            _s[np as usize - (n - delta)] = 0;
        }
    }
    return _s;
}

#[allow(non_snake_case)]
pub fn b_seq_piece(
    n: usize,
    delta: usize,
    m: usize,
    r: usize,
    qout: usize,
    s: &Vec<usize>,
) -> Vec<usize> {
    //index values that saves the index if a value in s gets sieved out (set 0)
    let mut index = Vec::new();
    let m0 = m + r;
    let mp = (m + 2 * r) as i64;
    let alpha0 = Rat {
        num: n % m0,
        den: m0,
    };
    let alpha1 = Rat {
        num: m0 * m0 - n % (m0 * m0),
        den: m0 * m0,
    };

    let eta = (Rational::from(delta) / Rational::from(m))
        * (Rational::from_f64(1.0).unwrap()
            + Rational::from_f64(1.0).unwrap() / Rational::from(qout));
    let (ainv, q) = dioph_appr(alpha1, 2 * r);
    let etaq = eta * Rational::from(q);
    let k = etaq.floor().to_f64() as usize;
    let c = (alpha0.num * q as usize + m0 / 2) / m0;
    let cainv = (ainv * c as i64) % q;
    let mut r0 = -cainv;
    for _j in 0..=k + 1 {
        if r0 <= -q {
            r0 += q;
        }
        let mut miter = m0 as i64 + r0;
        while miter >= m as i64 {
            if miter % 2 == 1 {
                let np = miter as usize * ((n + delta) / miter as usize);
                if np >= (n - delta)
                    && np <= (n + delta)
                    && np as i64 > miter
                    && s[np as usize - (n - delta)] == 1
                {
                    index.push(np as usize - (n - delta));
                    //s[np as usize - (n - delta)] = 0;
                }
            }
            miter -= q;
        }
        miter = m0 as i64 + r0 + q;
        while miter <= mp {
            if miter % 2 == 1 {
                let np = miter as usize * ((n + delta) / miter as usize);
                if np >= (n - delta)
                    && np <= (n + delta)
                    && np as i64 > miter
                    && s[np as usize - (n - delta)] == 1
                {
                    index.push(np as usize - (n - delta));
                    //s[np as usize - (n - delta)] = 0;
                }
            }
            miter += q;
        }
        r0 -= ainv;
    }

    r0 = -cainv + ainv;
    for _j in 0..k + 1 {
        if r0 > 0 {
            r0 -= q;
        }
        let mut miter = m0 as i64 + r0;
        while miter >= m as i64 {
            if miter % 2 == 1 {
                let np = miter as usize * ((n + delta) / miter as usize);
                if np >= (n - delta)
                    && np <= (n + delta)
                    && np as i64 > miter
                    && s[np as usize - (n - delta)] == 1
                {
                    index.push(np as usize - (n - delta));
                    //s[np as usize - (n - delta)] = 0;
                }
            }
            miter -= q;
        }
        miter = m0 as i64 + r0 + q;
        while miter <= mp {
            if miter % 2 == 1 {
                let np = miter as usize * ((n + delta) / miter as usize);
                if np >= (n - delta)
                    && np <= (n + delta)
                    && np as i64 > miter
                    && s[np as usize - (n - delta)] == 1
                {
                    index.push(np as usize - (n - delta));
                    //s[np as usize - (n - delta)] = 0;
                }
            }
            miter += q;
        }
        r0 += ainv;
    }
    return index;
}
