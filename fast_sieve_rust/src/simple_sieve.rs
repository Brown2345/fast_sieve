use rug::Integer;

#[allow(non_snake_case)]
pub fn simple_sieve(N: usize) -> Vec<i32>
// assumes an array of length >= N+1 is allocated at p
// ensures: for 1<=n<=N: p[n]=1 if n is prime, p[n]=0 otherwise
{
    let mut p = vec![0; N + 1];
    if N < 1 {
        return p;
    }
    if N < 2 {
        return p;
    }
    p[2] = 1;
    if N < 3 {
        return p;
    }

    let mut n = 3;

    while n <= N - 1 {
        p[n] = 1;
        n += 1;
        p[n] = 0;
        n += 1;
    }
    if n == N {
        p[n] = 1;
    }

    let mut m = 3;
    n = m * m;
    while n <= N {
        while n <= N {
            p[n] = 0;
            n += 2 * m;
        }
        m += 2;
        n = m * m;
    }
    return p;
}
//simple_sieve test here

#[allow(non_snake_case)]
pub fn simple_seg_sieve(n: usize, delta: usize, M: usize) -> Vec<usize>
//s has size D+1,
/* ensures, for 0<=j<=D, that
s[j] = 1   if n+j is coprime to all m<=M,
s[j] = 0   otherwise */
{
    let mut s = vec![1; delta + 1];

    //M<=1 means no prime number to check to
    if M <= 1 {
        for i in 0..=delta {
            s[i] = 1;
        }
        return s;
    }

    //we first set s[j]=1 for n+j odd, else zero
    let bn = n % 2;
    let bneq = (bn + 1) % 2;
    let mut j = 0;
    while j <= delta - 1 {
        s[j] = bn;
        j += 1;
        s[j] = bneq;
        j += 1;
    }
    if j == delta {
        s[j] = bn;
    }

    //special values depending on interval start being 0
    if n == 0 {
        s[0] = 0;
        if delta >= 1 {
            s[1] = 0;
        }
        if delta >= 2 {
            s[2] = 1;
        }
    }

    //special values depending on interval start being 1
    if n == 1 {
        s[0] = 0;
        if delta >= 1 {
            s[1] = 1;
        }
    }

    let smallprime = simple_sieve(M);
    let mut np;
    for m in (3..=M).step_by(2) {
        if smallprime[m] == 1 {
            np = m * (n as f64 / m as f64).ceil() as usize; //smallest multiple >=n of m
            if np <= m {
                np = 2 * m;
            }
            if np % 2 == 0 {
                np += m;
            }
            while np <= n + delta {
                s[np - n] = 0;
                np += 2 * m;
            }
        }
    }
    return s;
}

//test simple_seg_sieve

#[allow(non_snake_case)]
pub fn sub_seg_sieve(n: usize, delta: usize, M: usize) -> Vec<usize>
/* ensures, for 0<=j<=D, that
    s[j] = 1   if n+j is coprime to all m<=M,
    s[j] = 0   otherwise */
{
    let mut s = vec![0; delta + 1];

    //M<=1 means no prime number to check to
    if M <= 1 {
        for i in 0..=delta {
            s[i] = 1;
        }
        return s;
    }

    //we first set s[j]=1 for n+j odd, else zero
    let bn = n % 2;
    let bneq = (bn + 1) % 2;
    let mut j = 0;
    while j <= delta - 1 {
        s[j] = bn;
        j = j + 1;
        s[j] = bneq;
        j = j + 1;
    }
    if j == delta {
        s[j] = bn;
    }

    //special values depending on interval start being 0
    if n == 0 {
        if delta >= 1 {
            s[1] = 0;
        }
        if delta >= 2 {
            s[2] = 1;
        }
    }

    //special values depending on interval start being 1
    if n == 1 {
        s[0] = 0;
        if delta >= 1 {
            s[1] = 1;
        }
    }

    let deltap = Integer::from(M + 1).sqrt().to_usize_wrapping();
    let mut sqrt_md: usize;

    for m in (1..=M).step_by(deltap + 1) {
        sqrt_md = Integer::from(m + deltap).sqrt().to_usize_wrapping();
        let p = simple_seg_sieve(m, deltap, sqrt_md);

        let mut prime;
        if m % 2 == 1 {
            prime = m;
        } else {
            prime = m + 1;
        }
        while (prime <= m + deltap) && (prime <= M) {
            if p[prime - m] == 1 {
                let mut np = prime * (n as f64 / prime as f64).ceil() as usize; //smallest multiple >=n of m

                if np <= prime {
                    np = 2 * prime;
                }
                if np % 2 == 0 {
                    np += prime;
                }

                while np <= n + delta {
                    s[np - n] = 0;
                    np += 2 * prime;
                }
            }
            prime += 2;
        }
    }
    return s;
}

#[allow(non_snake_case)]
pub fn seg_sieve(n: usize, delta: usize)
//main function to call and to time
{
    let M = Integer::from(n + delta).sqrt().to_usize_wrapping();
    let _s = simple_seg_sieve(n, delta, M);
    println!("{:?}", _s);
}
