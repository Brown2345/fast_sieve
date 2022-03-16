pub struct Rat {
    pub num: i64,
    pub den: i64,
}

pub fn diophappr(mut alpha: Rat, Q: i64, a: &mut i64, ainv: &mut i64, qout: &mut i64)
/* precondition: alpha.num>=0, alpha.den>=0 */
/* constructs approximation a/q, q<=Q, to alpha */
/* sets ainv to a^{-1} mod q */
{
    let mut b: i64;
    let mut p: i64;
    let mut q: i64;
    let mut pmin: i64;
    let mut qmin: i64;
    let mut pplus: i64;
    let mut qplus: i64;
    let mut nummodb: i64;
    let mut flip: i64;
    let mut s: i32;

    //initial values
    b = alpha.num / alpha.den;
    p = b; pmin = 1; s = 1;
    q = 1; qmin = 0;

    while q <= Q {
        nummodb = alpha.num % alpha.den;
        //if alpha is a integer
        if nummodb == 0 {
            flip = if s==1 {-1*qmin} else {qmin};
            *a = p;
            ainv = if flip >= 0 {&mutflip % q} else {(q + flip) % q};
            *qout = q;
            return
        }
        alpha.num = alpha.den;
        alpha.den = nummodb;
        b = alpha.num / alpha.den;
        pplus = b * p + pmin;
        qplus = b * q + qmin;
        pmin = p;
        qmin = q;
        p = pplus;
        q = qplus;
        s = -s;
    }
    flip = if s == 1 {q} else {!q};
    *a = pmin;
    ainv = if flip >= 0 {flip % qmin} else {(qmin + flip) % qmin};
    *qout = qmin;
}
 fn main (){}
