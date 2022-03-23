pub fn small_sieve(miter: i64, n: usize, delta: usize) -> Vec<usize> {
    //sieves on a small interval to parallelize code
    let mut _s = vec![1; 2 * delta + 1];
    if miter % 2 == 1 {
        let np = ((n + delta) as f64 / miter as f64).floor() as i64 * miter;
        if np >= (n - delta) as i64 && np <= (n + delta) as i64 && np > miter {
            _s[np as usize - (n - delta)] = 0;
        }
    }
    return _s;
}
