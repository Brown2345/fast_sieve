pub struct Rat {
    pub num: usize,
    pub den: usize,
}

#[allow(non_snake_case)]
pub fn dioph_appr(mut alpha: Rat, Q: usize) -> (i64, i64)
/* precondition: alpha.num>=0, alpha.den>=0 */
/* constructs approximation a/q, q<=Q, to alpha */
/* sets ainv to a^{-1} mod q */
{
    let mut b = alpha.num / alpha.den;
    let mut p = b;
    let mut q: i64 = 1;
    let mut pmin = 1;
    let mut qmin = 0;
    let mut s = 1 as i64;
    let ainv: i64;
    let qout: i64;

    while q <= Q as i64 {
        let nummodb = alpha.num % alpha.den;
        if nummodb == 0 {
            if -s * qmin >= 0 {
                ainv = (-s * qmin) % q;
            } else {
                ainv = (q - 1) - (-s * qmin - 1) % q;
            }
            qout = q;
            return (ainv, qout);
        }
        alpha.num = alpha.den;
        alpha.den = nummodb;
        b = alpha.num / alpha.den;
        let pplus = b * p + pmin;
        let qplus = b as i64 * q + qmin;
        pmin = p;
        qmin = q;
        p = pplus;
        q = qplus;
        s = -s;
    }
    if s * q >= 0 {
        ainv = (s * q) % qmin;
    } else {
        ainv = (qmin - 1) - (-s * q - 1) % qmin;
    }
    qout = qmin;
    return (ainv, qout);
}

//test here
