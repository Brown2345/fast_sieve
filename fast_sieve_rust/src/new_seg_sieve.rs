use crate::dioph_appr::dioph_appr;
use crate::simple_sieve::sub_seg_sieve;

use crate::dioph_appr::Rat;

/// algorithm 1 from Helfgott, returns array with values 1 in index j if n+j is a coprime
/// sequential version
#[allow(non_snake_case)]
pub fn new_seg_sieve(n: usize, delta: usize, k: usize, quot: usize) -> Vec<usize>
// Algorithm 1 from Helfgott, takes input n, delta, K, and ..
{
    //return binary vector s with 0 and 1 depending of the index is a prime or not
    let mut s = vec![0; 2 * delta + 1];

    let sqrt_n_delta = (((n + delta) as f64).sqrt()).floor() as usize;
    s = sub_seg_sieve(n - delta, 2 * delta, k * delta);

    let mut r: usize = 0;
    for m in ((k * delta + 1)..=sqrt_n_delta).step_by(2 * r + 1) {
        let m2kap = m * m * delta / (quot * n);
        r = ((m2kap as f64).sqrt()).floor() as usize;

        let m0 = m + r;
        let mp = m + 2 * r;
        let alpha0 = Rat {
            num: n % m0,
            den: m0,
        };
        let alpha1 = Rat {
            num: m0 * m0 - n % (m0 * m0),
            den: m0 * m0,
        };
        let eta = (delta as f64 / m as f64) * (1.0 + 1.0 / quot as f64);
        let a = 0;
        let ainv = 0;
        let q = 0;
        dioph_appr(alpha1, 2 * r, a, ainv, q);
        let k = (eta * q as f64).floor() as usize;
        let c = (alpha0.num * q as usize + m0 / 2) / m0;
        let cainv = (ainv * c as i64) % q as i64;
        let mut r0 = -cainv;
        for _j in 0..=k + 1 {
            if r0 <= -(q as i64) {
                r0 = r0 + q as i64;
            }
            let mut miter = m0 + r0 as usize;
            while miter >= m {
                if miter % 2 == 0 {
                    let np = ((n + delta) as f64 / miter as f64).floor() as usize * miter;
                    if np >= n - delta && np <= n + delta && np > miter {
                        s[np - (n - delta)] = 0;
                    }
                }
                miter = miter - q as usize;
            }
            miter = m0 + r0 as usize + q as usize;
            while miter <= mp {
                if miter % 2 == 0 {
                    let np = ((n + delta) as f64 / miter as f64).floor() as usize * miter;
                    if np >= n - delta && np <= n + delta && np > miter {
                        s[np - (n - delta)] = 0;
                    }
                }
                miter = miter + q as usize;
            }
            r0 = r0 - ainv;
        }

        r0 = -cainv + ainv;
        for _j in 0..k + 1 {
            if r0 > 0 {
                r0 = r0 - q;
            }
            let mut miter = m0 + r0 as usize;
            while miter >= m {
                if miter % 2 == 0 {
                    let np = ((n + delta) as f64 / miter as f64).floor() as usize * miter;
                    if np >= n - delta && np <= n + delta && np > miter {
                        s[np - (n - delta)] = 0;
                    }
                }
                miter = miter - q as usize;
            }
            miter = m0 + r0 as usize + q as usize;
            while miter <= mp {
                if miter % 2 == 0 {
                    let np = ((n + delta) as f64 / miter as f64).floor() as usize * miter;
                    if np >= n - delta && np <= n + delta && np > miter {
                        s[np - (n - delta)] = 0;
                    }
                }
                miter = miter + q as usize;
            }
            r0 = r0 + ainv;
        }
    }
    return s;
}
