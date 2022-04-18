#![allow(unused_imports)]
use fast_sieve_rust::par_pieces::b_seq_piece;
use fast_sieve_rust::par_sieve_atomic::par_sieve_atomic;
use fast_sieve_rust::parallel;
use fast_sieve_rust::parallel::parallel_sieve;
use fast_sieve_rust::s_par_sieve;
use fast_sieve_rust::sieve_loop;
use fast_sieve_rust::sieve_s_iter;
use fast_sieve_rust::simple_seg_sieve;
use rayon::prelude::*;
use std::time::Instant;

fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(6)
        .build_global()
        .unwrap();

    let n: usize = 5000000000000000000;
    let delta = 40000000;
    let time = Instant::now();
    let _s = sieve_loop(n, delta, 8, 4);
    let elapsed_time = time.elapsed();
    println!("Sequential is = {}", elapsed_time.as_millis());
    let time2 = Instant::now();
    let _s = parallel_sieve(n, delta, 8, 4);
    let elapsed_time2 = time2.elapsed();
    println!("Parallel is = {}", elapsed_time2.as_millis());
}
