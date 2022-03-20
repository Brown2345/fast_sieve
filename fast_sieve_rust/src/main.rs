use fast_sieve_rust::seg_sieve;
use fast_sieve_rust::simple_seg_sieve;
use fast_sieve_rust::sub_seg_sieve;

fn main() {
    let s = seg_sieve(1000, 100);
    println!("Seg is {:?}", s);
    for i in 0..=100 {
        if s[i] == 1 {
            println!("Eine Primzahl ist {}", 1000 + i);
        }
    }
}
