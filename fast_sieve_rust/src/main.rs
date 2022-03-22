#![allow(unused_imports)]
use fast_sieve_rust::new_seg_sieve;
use fast_sieve_rust::seg_sieve;
use fast_sieve_rust::simple_seg_sieve;
use fast_sieve_rust::simple_sieve;
use fast_sieve_rust::sub_seg_sieve;

fn main() {
    let n = 100;
    let delta = 50;
    let s = new_seg_sieve(n, delta, 8, 4);
    println!("Seg is {:?}", s);
    for i in 0..=2 * delta {
        if s[i] == 1 {
            println!("Eine Primzahl ist {}", n - delta + i);
        }
    }
}
