#![allow(unused_imports)]
use crate::dioph_appr::{dioph_appr, Rat};
use crate::par_pieces::b_seq_piece;
use crate::par_pieces_atomic::{b_seq_piece_atomic, sub_seg_sieve_atomic};
use crate::simple_sieve::sub_seg_sieve;
use core::f64;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::prelude::*;
use rug::{Integer, Rational};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

#[derive(Debug, Clone)]
struct SimpleBitVec(Box<[u32]>);
impl SimpleBitVec {
    fn new_set(n: usize) -> Self {
        let len = (n + 31) / 32;
        SimpleBitVec(
            (0..len)
                .map(|_| 0xFFFFFFFF)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        )
    }
    fn unset(&mut self, i: usize) {
        let (idx, bit) = (i / 32, i % 32);
        self.0[idx] &= !(1 << bit);
    }
    fn intersect_mut(&mut self, other: &Self) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a &= *b;
        }
    }
    fn get(&self, i: usize) -> u8 {
        let (idx, bit) = (i / 32, i % 32);
        let value = self.0[idx] & (1 << bit);
        if value != 0 {
            1
        } else {
            0
        }
    }
}

pub fn parallel() {
    let modulus: u32 = 7000000;
    let count: u32 = 100000000;

    let mut rng = StdRng::from_entropy();
    println!("a");
    let sim_input = (0..count)
        .map(move |_| rng.gen_range(0..100 * modulus))
        .collect::<Vec<_>>();
    println!("a");
    println!("{:?}", sim_input);
    let mk_new_bitvec = move || SimpleBitVec::new_set(modulus as usize);
    let result = sim_input
        .par_iter()
        .fold(mk_new_bitvec, move |mut acc, &x| {
            // This is where you'd put your update logic; this
            // example code just randomly unsets a bit 1% of the time
            if x < modulus {
                acc.unset(x as usize);
            }
            acc
        })
        .reduce(mk_new_bitvec, |mut acc, x| {
            acc.intersect_mut(&x);
            acc
        });

    println!(
        "{:?}",
        (0..modulus as usize)
            .map(|i| result.get(i))
            .collect::<Vec<_>>()
    );
    //here return s
}

/// algorithm 1 from Helfgott, returns array with values 1 in index j if n+j is a coprime
/// sequential version
pub fn par_sieve_disc(n: usize, delta: usize, k: usize, qout: usize) -> Vec<usize>
// Algorithm 1 from Helfgott, takes input n, delta, K, and ..
{
    //return binary vector s with 0 and 1 depending of the index is a prime or not

    let nplusd = Integer::from(n + delta);
    let sqrt_n_delta = nplusd.sqrt().to_usize_wrapping();
    let s = sub_seg_sieve(n - delta, 2 * delta, k * delta);
    //length and important values for the big loop
    let mut _count = 0;
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
        _count += 1;
        m += 2 * r + 1;
    }

    parallel();
    return s;
}
