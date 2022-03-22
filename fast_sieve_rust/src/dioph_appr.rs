pub struct Rat {
    pub num: usize,
    pub den: usize,
}

#[allow(non_snake_case)]
pub fn dioph_appr(mut alpha: Rat, Q: usize, mut _a: usize, mut _ainv: i64, mut _qout: i64)
/* precondition: alpha.num>=0, alpha.den>=0 */
/* constructs approximation a/q, q<=Q, to alpha */
/* sets ainv to a^{-1} mod q */
{
    let mut b = alpha.num / alpha.den;
    let mut p = b;
    let mut q = 1;
    let mut pmin = 1;
    let mut qmin = 0;
    let mut s = 1 as i64;

    while q <= Q {
        let nummodb = alpha.num % alpha.den;
        if nummodb == 0 {
            _a = p;
            _ainv = (-s * qmin as i64) % q as i64;
            _qout = q as i64;
            return;
        }
        alpha.num = alpha.den;
        alpha.den = nummodb;
        b = alpha.num / alpha.den;
        let pplus = b * p + pmin;
        let qplus = b * q + qmin;
        pmin = p;
        qmin = q;
        p = pplus;
        q = qplus;
        s = -s;
    }
    _a = pmin;
    _ainv = (s * q as i64) % qmin as i64;
    _qout = qmin as i64;
}

//test here
