#![allow(unused_imports)]
use crate::dioph_appr::{dioph_appr, Rat};
use crate::par_pieces::b_seq_piece;
use crate::par_pieces_atomic::{b_seq_piece_atomic, sub_seg_sieve_atomic};
use crate::simple_sieve::sub_seg_sieve;
use core::f64;
use rayon::prelude::*;
use rug::{Integer, Rational};

/// algorithm 1 from Helfgott, returns array with values 1 in index j if n+j is a coprime

pub fn parallel_sieve(n: usize, delta: usize, k: usize, qout: usize) -> Vec<usize>
// Algorithm 1 from Helfgott, takes input n, delta, K, and ..
{
    //return binary vector s with 0 and 1 depending of the index is a prime or not

    let nplusd = Integer::from(n + delta);
    let sqrt_n_delta = nplusd.sqrt().to_usize_wrapping();
    let s = sub_seg_sieve(n - delta, 2 * delta, k * delta);
    //length and important values for the big loop
    let mut count = 0;
    let mut r_val = Vec::new();
    let mut m_val = Vec::new();
    //counting the length of the big loop
    let mut r;
    let mut m = k * delta + 1;
    while m <= sqrt_n_delta {
        let m2kap = Integer::from(m) * Integer::from(m) * Integer::from(delta)
            / (Integer::from(qout) * Integer::from(n)); //not sure if correct
        r = m2kap.sqrt().to_usize_wrapping();
        r_val.push(r);
        m_val.push(m);
        count += 1;
        m += 2 * r + 1;
    }
    let num_threads = 6;
    let h = count / num_threads;
    let borders = vec![
        vec![0, h],
        vec![h + 1, 2 * h],
        vec![2 * h + 1, 3 * h],
        vec![3 * h + 1, 4 * h],
        vec![4 * h + 1, 5 * h],
        vec![5 * h + 1, count - 1],
    ];
    let mut results = vec![s.clone(); num_threads];
    (0..num_threads)
        .into_par_iter()
        .zip(results.par_iter_mut())
        .for_each(|(i, interval)| {
            let left = borders[i][0];
            let right = borders[i][1];
            *interval = sub_interval(
                n,
                delta,
                m_val.clone(),
                r_val.clone(),
                qout,
                s.clone(),
                left,
                right,
            )
        });
    let sum: Vec<usize> = results
        .into_iter()
        .reduce(|a, b| a.iter().zip(&b).map(|(x, y)| x * y).collect())
        .unwrap();
    return sum;
}

pub fn sub_interval(
    n: usize,
    delta: usize,
    m: Vec<usize>,
    r: Vec<usize>,
    qout: usize,
    mut s: Vec<usize>,
    left: usize,
    right: usize,
) -> Vec<usize> {
    for iter in left..=right {
        let m0 = m[iter] + r[iter];
        let mp = (m[iter] + 2 * r[iter]) as i64;
        let alpha0 = Rat {
            num: n % m0,
            den: m0,
        };
        let alpha1 = Rat {
            num: m0 * m0 - n % (m0 * m0),
            den: m0 * m0,
        };
        let eta = (Rational::from(delta) / Rational::from(m[iter]))
            * (Rational::from_f64(1.0).unwrap()
                + Rational::from_f64(1.0).unwrap() / Rational::from(qout));
        let (ainv, q) = dioph_appr(alpha1, 2 * r[iter]);
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
            while miter >= m[iter] as i64 {
                if miter % 2 == 1 {
                    let np = miter as usize * ((n + delta) / miter as usize);
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        s[np as usize - (n - delta)] = 0;
                    }
                }
                miter -= q;
            }
            miter = m0 as i64 + r0 + q;
            while miter <= mp {
                if miter % 2 == 1 {
                    let np = miter as usize * ((n + delta) / miter as usize);
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        s[np as usize - (n - delta)] = 0;
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
            while miter >= m[iter] as i64 {
                if miter % 2 == 1 {
                    let np = miter as usize * ((n + delta) / miter as usize);
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        s[np as usize - (n - delta)] = 0;
                    }
                }
                miter -= q;
            }
            miter = m0 as i64 + r0 + q;
            while miter <= mp {
                if miter % 2 == 1 {
                    let np = miter as usize * ((n + delta) / miter as usize);
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        s[np as usize - (n - delta)] = 0;
                    }
                }
                miter += q;
            }
            r0 += ainv;
        }
    }
    return s;
}
