/*
#![allow(unused_imports)]
use crate::dioph_appr::{dioph_appr, Rat};
use crate::par_pieces_ptr::par_ptr;
use crate::simple_sieve::sub_seg_sieve;
use core::f64;
use rayon::prelude::*;
use rug::{Integer, Rational};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

/// algorithm 1 from Helfgott, returns array with values 1 in index j if n+j is a coprime
/// sequential version
pub fn par_sieve_ptr(n: usize, delta: usize, k: usize, qout: usize) -> Vec<usize>
// Algorithm 1 from Helfgott, takes input n, delta, K, and ..
{
    //return binary vector s with 0 and 1 depending of the index is a prime or not

    let nplusd = Integer::from(n + delta);
    let sqrt_n_delta = nplusd.sqrt().to_usize_wrapping();
    let s = sub_seg_sieve(n - delta, 2 * delta, k * delta);

    let (s, len, cap) = s.into_raw_parts();

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
    iterator
        .into_par_iter()
        .for_each(|i: usize| par_ptr(n, delta, m_val[i], r_val[i], qout, s, len, cap));
    let result = unsafe { Vec::from_raw_parts(s, len, cap) };
    return result;
}
*/
