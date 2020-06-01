use ndarray::prelude::*;

#[allow(dead_code)]
pub(super) fn sproot(
    t: ArrayView1<f64>,
    n: usize,
    c: ArrayView1<f64>,
    mut zero: ArrayViewMut1<f64>,
    mest: usize,
    m: &mut usize,
    ier: &mut usize,
) {
    let mut y = Array1::<f64>::zeros(3);
    let n4 = n - 4;

    let mut a0 = 0.0;
    let mut a1 = 0.0;
    let mut a2 = 0.0;
    let mut a3 = 0.0;
    let mut b0 = 0.0;
    let mut b1 = 0.0;
    let mut c1 = 0.0;
    let mut c2 = 0.0;

    *ier = 10;
    if n < 8 {
        return;
    }
    let mut j = n;
    for i in 1..4 {
        if t[i - 1] > t[i + 1 - 1] || t[j - 1] < t[j - 1 - 1] {
            return;
        }
        j = j - 1;
    }
    for i in 4..(n4 + 1) {
        if t[i - 1] >= t[i + 1 - 1] {
            return;
        }
    }
    // the problem considered reduces to finding the zeros of the cubic
    // polynomials pl(x) which define the cubic spline in each knot
    // interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
    // the condition that it belongs to the knot interval.
    // the cubic polynomial pl(x) is determined by computing s(t(l)),
    // s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
    // s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
    // splines and their derivatives, the value of s(t(l)) and s'(t(l))
    // is already known from the foregoing knot interval.
    *ier = 0;
    // evaluate some constants for the first knot interval
    let mut h1 = t[4 - 1] - t[3 - 1];
    let mut h2 = t[5 - 1] - t[4 - 1];
    let mut t1 = t[4 - 1] - t[2 - 1];
    let mut t2 = t[5 - 1] - t[3 - 1];
    let mut t3 = t[6 - 1] - t[4 - 1];
    let mut t4 = t[5 - 1] - t[2 - 1];
    let mut t5 = t[6 - 1] - t[3 - 1];
    // calculate a0 = s(t(4)) and ah = s'(t(4)).
    c1 = c[1 - 1];
    c2 = c[2 - 1];
    let mut c3 = c[3 - 1];
    let mut c4 = (c2 - c1) / t4;
    let mut c5 = (c3 - c2) / t5;
    let mut d4 = (h2 * c1 + t1 * c2) / t4;
    let mut d5 = (t3 * c2 + h1 * c3) / t5;
    a0 = (h2 * d4 + h1 * d5) / t2;
    let mut ah = 3.0 * (h2 * c4 + h1 * c5) / t2;
    let mut z1 = true;
    if ah < 0.0 {
        z1 = false;
    }
    let mut nz1 = !z1;
    *m = 0;
    // main loop for the different knot intervals.
    for l in 4..(n4 + 1) {
        // evaluate some constants for the knot interval t(l) <= x <= t(l+1).
        h1 = h2;
        h2 = t[l + 2 - 1] - t[l + 1 - 1];
        t1 = t2;
        t2 = t3;
        t3 = t[l + 3 - 1] - t[l + 1 - 1];
        t4 = t5;
        t5 = t[l + 3 - 1] - t[l - 1];
        // find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
        c1 = c2;
        c2 = c3;
        c3 = c[l - 1];
        // c4 = c5;
        // c5 = (c3 - c2) / t5;
        d4 = (h2 * c1 + t1 * c2) / t4;
        d5 = (h1 * c3 + t3 * c2) / t5;
        let b0 = (h2 * d4 + h1 * d5) / t2;
        // let bh = 3.0 * (h2 * c4 + h1 * c5) / t2;
        // test whether pl(x) could have a zero in the range
        // t(l) <= x <= t(l+1).
        let mut z3 = true;
        if b1 < 0.00 {
            z3 = false;
        }
        let nz3 = !z3;
        if a0 * b0 > 0.0 {
            let mut z0 = true;
            if a0 < 0.0 {
                z0 = false;
            }
            let nz0 = !z0;
            let mut z2 = true;
            if a2 < 0.0 {
                z2 = false;
            }
            let nz2 = !z2;
            let mut z4 = true;
            if 3.0 * a3 + a2 < 0.0 {
                z4 = false;
            }
            let nz4 = !z4;
            if !((z0 && ((nz1 && (z3 || (z2 && nz4))) || (nz2 && z3 && z4)))
                || (nz0 && ((z1 && (nz3 || (nz2 && z4))) || (z2 && nz3 && nz4))))
            {
                a0 = b0;
                // ah = bh;
                z1 = z3;
                nz1 = nz3;
                continue;
            }
        }
        // find the zeros of ql(y).
        j = super::fpcuro::fpcuro(a3, a2, a1, a0, y.view_mut());
        if j == 0 {
            a0 = b0;
            // ah = bh;
            z1 = z3;
            nz1 = nz3;
            continue;
        }
        // find which zeros of pl(x) are zeros of s(x).
        for i in 1..(j + 1) {
            if y[i - 1] < 0.0 || y[i - 1] > 1.0 {
                continue;
            }
            //  test whether the number of zeros of s(x) exceeds mest.
            if *m >= mest {
                *ier = 1;
                return;
            }
            *m += 1;
            zero[*m - 1] = t[l - 1] + h1 * y[i - 1];
        }
        a0 = b0;
        // ah = bh;
        z1 = z3;
        nz1 = nz3;
    }
    // the zeros of s(x) are arranged in increasing order.
    if *m < 2 {
        return;
    }
    for i in 2..(*m + 1) {
        j = i;
        loop {
            let j1 = j - 1;
            if j1 == 0 {
                break;
            }
            if zero[j - 1] >= zero[j1 - 1] {
                break;
            }
            let zz = zero[j - 1];
            zero[j - 1] = zero[j1 - 1];
            zero[j1 - 1] = zz;
            j = j1;
        }
    }
    j = *m;
    *m = 1;
    for i in 2..(j + 1) {
        if zero[i - 1] == zero[*m - 1] {
            continue;
        }
        *m += 1;
        zero[*m - 1] = zero[i - 1];
    }
}
