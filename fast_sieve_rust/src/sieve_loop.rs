#![allow(unused_imports)]
use crate::dioph_appr::{dioph_appr, Rat};
use crate::simple_sieve::{simple_seg_sieve, sub_seg_sieve};

use core::f64;
use rug::{Integer, Rational};
//pub struct Mpz = mpz::new();
/// algorithm 1 from Helfgott, returns array with values 1 in index j if n+j is a coprime
/// sequential version
#[allow(non_snake_case)]
pub fn sieve_loop(n: usize, delta: usize, k: usize, quot: usize) -> Vec<usize>
// Algorithm 1 from Helfgott, takes input n, delta, K, and ..
{
    //return binary vector s with 0 and 1 depending of the index is a prime or not
    let mut _s = vec![0; 2 * delta + 1];

    let npD = Integer::from(n + delta);
    let sqrt_n_delta = npD.sqrt().to_usize_wrapping();
    _s = sub_seg_sieve(n - delta, 2 * delta, k * delta);

    let mut r;
    let mut m = k * delta + 1;
    while m <= sqrt_n_delta {
        let m2kap = Integer::from(m) * Integer::from(m) * Integer::from(delta)
            / (Integer::from(quot) * Integer::from(n)); //not sure if correct
        r = m2kap.sqrt().to_usize_wrapping();
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
                + Rational::from_f64(1.0).unwrap() / Rational::from(quot));
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
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        _s[np as usize - (n - delta)] = 0;
                    }
                }
                miter -= q;
            }
            miter = m0 as i64 + r0 + q;
            while miter <= mp {
                if miter % 2 == 1 {
                    let np = miter as usize * ((n + delta) / miter as usize);
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        _s[np as usize - (n - delta)] = 0;
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
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        _s[np as usize - (n - delta)] = 0;
                    }
                }
                miter -= q;
            }
            miter = m0 as i64 + r0 + q;
            while miter <= mp {
                if miter % 2 == 1 {
                    let np = miter as usize * ((n + delta) / miter as usize);
                    if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                        _s[np as usize - (n - delta)] = 0;
                    }
                }
                miter += q;
            }
            r0 += ainv;
        }
        m += 2 * r + 1;
    }
    return _s;
}
