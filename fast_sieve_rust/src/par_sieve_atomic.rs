#![allow(unused_imports)]
use crate::dioph_appr::{dioph_appr, Rat};
use crate::par_pieces::b_seq_piece;
use crate::par_pieces_atomic::{b_seq_piece_atomic, sub_seg_sieve_atomic};
use crate::simple_sieve::sub_seg_sieve;
use core::f64;
use rayon::prelude::*;
use rug::{Integer, Rational};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

/// algorithm 1 from Helfgott, returns array with values 1 in index j if n+j is a coprime
/// sequential version
pub fn par_sieve_atomic(n: usize, delta: usize, k: usize, qout: usize) -> Vec<AtomicUsize>
// Algorithm 1 from Helfgott, takes input n, delta, K, and ..
{
    //return binary vector s with 0 and 1 depending of the index is a prime or not

    let nplusd = Integer::from(n + delta);
    let sqrt_n_delta = nplusd.sqrt().to_usize_wrapping();
    let s = sub_seg_sieve_atomic(n - delta, 2 * delta, k * delta);
    let s_arc = Arc::new(s);
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
    let iterator = 0..count;
    let mut threads = Vec::with_capacity(count);
    for i in iterator {
        let s_arc = Arc::clone(&s_arc);
        let m_val_i = m_val[i];
        let r_val_i = r_val[i];
        threads.push(thread::spawn(move || {
            b_seq_piece_atomic(n, delta, m_val_i, r_val_i, qout, s_arc)
        }));
    }
    for thread in threads {
        thread.join().unwrap();
    }
    return Arc::try_unwrap(s_arc).unwrap();
}
