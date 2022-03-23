#![allow(unused_imports)]
use fast_sieve_rust::new_seg_sieve;
use fast_sieve_rust::parallel_sieve;
use fast_sieve_rust::seg_sieve;
use fast_sieve_rust::simple_seg_sieve;
use fast_sieve_rust::simple_sieve;
use fast_sieve_rust::small_sieve;
use fast_sieve_rust::sub_seg_sieve;
use rayon::prelude::*;

fn main() {
    /*
    let n = 10;
    let delta = 5;
    let deltak = 40;
    let x = 5;
    let y = 6;
    let z = 1;
    let iterator = (x..=y).step_by(z);
    let len = iterator.clone().count();
    let mut returns = vec![vec![0; deltak]; len];
    //let startvec = vec![1; iterator.clone().count()];
    iterator
        .into_iter()
        .zip(returns.iter_mut())
        .par_bridge()
        .for_each(|(i, r)| *r = small_sieve(i, n, delta));
    let _result = returns.into_par_iter().reduce(
        || vec![1; len],
        |a, b| a.iter().zip(&b).map(|(x, y)| x * y).collect(),
    );
    */
    let n = 500;
    let delta = 100;
    let s = parallel_sieve(n, delta, 8, 4);
    println!("Seg is {:?}", s);
    for i in 0..=2 * delta {
        if s[i] == 1 {
            println!("Eine Primzahl ist {}", n - delta + i);
        }
    }
}
