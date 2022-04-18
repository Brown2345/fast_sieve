#![allow(unused_imports)]
use crate::dioph_appr::{dioph_appr, Rat};
use crate::simple_sieve::simple_seg_sieve;
use rug::{Integer, Rational};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

#[allow(non_snake_case)]
pub fn sub_seg_sieve_atomic(n: usize, delta: usize, M: usize) -> Vec<AtomicUsize>
/* ensures, for 0<=j<=D, that
    s[j] = 1   if n+j is coprime to all m<=M,
    s[j] = 0   otherwise */
{
    let length = 2 * delta + 1;
    let s: Vec<AtomicUsize> = std::iter::repeat_with(|| AtomicUsize::new(0))
        .take(length)
        .collect();
    //M<=1 means no prime number to check to
    if M <= 1 {
        for i in 0..=delta {
            s[i].store(1, Ordering::Relaxed);
        }
        return s;
    }

    //we first set s[j]=1 for n+j odd, else zero
    let bn = n % 2;
    let bneq = (bn + 1) % 2;
    let mut j = 0;
    while j <= delta - 1 {
        s[j].store(bn, Ordering::Relaxed);
        j = j + 1;
        s[j].store(bneq, Ordering::Relaxed);
        j = j + 1;
    }
    if j == delta {
        s[j].store(bn, Ordering::Relaxed);
    }

    //special values depending on interval start being 0
    if n == 0 {
        if delta >= 1 {
            s[1].store(0, Ordering::Relaxed);
        }
        if delta >= 2 {
            s[2].store(1, Ordering::Relaxed);
        }
    }

    //special values depending on interval start being 1
    if n == 1 {
        s[0].store(0, Ordering::Relaxed);
        if delta >= 1 {
            s[1].store(1, Ordering::Relaxed);
        }
    }

    let deltap = Integer::from(M + 1).sqrt().to_usize_wrapping();
    let mut sqrt_md: usize;

    for m in (1..=M).step_by(deltap + 1) {
        sqrt_md = Integer::from(m + deltap).sqrt().to_usize_wrapping();
        let p = simple_seg_sieve(m, deltap, sqrt_md);

        let mut prime;
        if m % 2 == 1 {
            prime = m;
        } else {
            prime = m + 1;
        }
        while (prime <= m + deltap) && (prime <= M) {
            if p[prime - m] == 1 {
                let mut np = prime * ((n + prime - 1) / prime); //smallest multiple >=n of m

                if np <= prime {
                    np = 2 * prime;
                }
                if np % 2 == 0 {
                    np += prime;
                }

                while np <= n + delta {
                    s[np - n].store(0, Ordering::Relaxed);
                    np += 2 * prime;
                }
            }
            prime += 2;
        }
    }
    return s;
}

#[allow(non_snake_case)]
pub fn b_seq_piece_atomic(
    n: usize,
    delta: usize,
    m: usize,
    r: usize,
    qout: usize,
    s: Arc<Vec<AtomicUsize>>,
) {
    //index values that saves the index if a value in s gets sieved out (set 0)
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
                if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                    s[np as usize - (n - delta)].store(0, Ordering::Relaxed);
                }
            }
            miter -= q;
        }
        miter = m0 as i64 + r0 + q;
        while miter <= mp {
            if miter % 2 == 1 {
                let np = miter as usize * ((n + delta) / miter as usize);
                if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                    s[np as usize - (n - delta)].store(0, Ordering::Relaxed);
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
                    s[np as usize - (n - delta)].store(0, Ordering::Relaxed);
                }
            }
            miter -= q;
        }
        miter = m0 as i64 + r0 + q;
        while miter <= mp {
            if miter % 2 == 1 {
                let np = miter as usize * ((n + delta) / miter as usize);
                if np >= (n - delta) && np <= (n + delta) && np as i64 > miter {
                    s[np as usize - (n - delta)].store(0, Ordering::Relaxed);
                }
            }
            miter += q;
        }
        r0 += ainv;
    }
}
