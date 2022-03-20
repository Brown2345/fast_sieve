pub fn new_seg_sieve(n: usize, delta: usize) -> Vec<usize>
// Algorithm 1 from Helfgott, takes input n, delta, K, and ..
{
    //return binary vector s with 0 and 1 depending of the index is a prime or not
    let mut s = vec![0; 2 * delta + 1];

    let sqrt_n_delta = (((n + delta) as f64).sqrt()).floor();
    s = sub_seg_sieve(...);
    return s;
}
