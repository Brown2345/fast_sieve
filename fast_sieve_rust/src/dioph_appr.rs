#[allow(non_snake_case)]
pub fn dioph_appr(mut alpha: f64, Q: i64, a: &mut i64, ainv: &mut i64, qout: &mut i64)
/* precondition: alpha.num>=0, alpha.den>=0 */
/* constructs approximation a/q, q<=Q, to alpha */
/* sets ainv to a^{-1} mod q */
{
    let mut b: i64 = alpha.floor() as i64;
    let mut p: i64 = b;
    let mut q: i64 = 1;
    let mut pmin: i64 = 1;
    let mut qmin: i64 = 0;
    let mut pplus: i64;
    let mut qplus: i64;
    let mut s: i64 = 1;

    while q <= Q  {
        if alpha == b as f64{
            *a = p;
            *ainv = (-s * qmin) % q;
            *qout = q;
            return
        }
        alpha = 1.0 / (alpha - b as f64);
        b = alpha.floor() as i64;
        pplus = b * p + pmin;
        qplus = b * q + qmin;
        pmin = p;
        qmin = q;
        p = pplus;
        q = qplus;
        s = -s;
    }
    *a = pmin;
    *ainv = (s * q) % qmin;
    *qout = qmin;
}

//test here
