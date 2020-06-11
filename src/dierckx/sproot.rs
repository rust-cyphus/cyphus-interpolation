use ndarray::prelude::*;

#[allow(dead_code)]
pub(super) fn sproot(
    t: ArrayView1<f64>,
    n: usize,
    c: ArrayView1<f64>,
    ier: &mut usize,
) -> Array1<f64> {
    let mut zero = Vec::<f64>::new();
    let mut y = Array1::<f64>::zeros(3);
    let n4 = n - 4;

    *ier = 10;
    if n < 8 {
        return Array::from(zero);
    }
    let mut j = n;
    for i in 1..=3 {
        if t[i - 1] > t[i] || t[j - 1] < t[j - 2] {
            return Array::from(zero);
        }
        j -= 1;
    }
    for i in 4..=n4 {
        if t[i - 1] >= t[i] {
            return Array::from(zero);
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
    let mut h1 = t[3] - t[2];
    let mut h2 = t[4] - t[3];

    let mut t1 = t[3] - t[1];
    let mut t2 = t[4] - t[2];
    let mut t3 = t[5] - t[3];
    let mut t4 = t[4] - t[1];
    let mut t5 = t[5] - t[2];

    // calculate a0 = s(t(4)) and ah = s'(t(4)).
    let mut c1 = c[0];
    let mut c2 = c[1];
    let mut c3 = c[2];
    let mut c4 = (c2 - c1) / t4;
    let mut c5 = (c3 - c2) / t5;

    let mut d4 = (h2 * c1 + t1 * c2) / t4;
    let mut d5 = (t3 * c2 + h1 * c3) / t5;

    let mut a0 = (h2 * d4 + h1 * d5) / t2;
    let mut ah = 3.0 * (h2 * c4 + h1 * c5) / t2;
    let mut z1 = ah >= 0.0;
    let mut nz1 = !z1;
    // main loop for the different knot intervals.
    for l in 4..=n4 {
        // evaluate some constants for the knot interval t(l) <= x <= t(l+1).
        h1 = h2;
        h2 = t[l + 1] - t[l];
        t1 = t2;
        t2 = t3;
        t3 = t[l + 2] - t[l];
        t4 = t5;
        t5 = t[l + 2] - t[l - 1];
        // find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
        c1 = c2;
        c2 = c3;
        c3 = c[l - 1];
        c4 = c5;
        c5 = (c3 - c2) / t5;
        d4 = (h2 * c1 + t1 * c2) / t4;
        d5 = (h1 * c3 + t3 * c2) / t5;
        let b0 = (h2 * d4 + h1 * d5) / t2;
        let bh = 3.0 * (h2 * c4 + h1 * c5) / t2;
        // calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
        // pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
        let a1 = ah * h1;
        let b1 = bh * h1;
        let a2 = 3.0 * (b0 - a0) - b1 - 2.0 * a1;
        let a3 = 2.0 * (a0 - b0) + b1 + a1;
        // test whether pl(x) could have a zero in the range
        // t(l) <= x <= t(l+1).
        let z3 = b1 >= 0.0;
        let nz3 = !z3;
        if a0 * b0 > 0.0 {
            let z0 = a0 >= 0.0;
            let nz0 = !z0;
            let z2 = a2 >= 0.0;
            let nz2 = !z2;
            let z4 = 3.0 * a3 + a2 >= 0.0;
            let nz4 = !z4;
            if !((z0 && ((nz1 && (z3 || (z2 && nz4))) || (nz2 && z3 && z4)))
                || (nz0 && ((z1 && (nz3 || (nz2 && z4))) || (z2 && nz3 && nz4))))
            {
                a0 = b0;
                ah = bh;
                z1 = z3;
                nz1 = nz3;
                continue;
            }
        }
        // find the zeros of ql(y).
        j = super::fpcuro::fpcuro(a3, a2, a1, a0, y.view_mut());
        if j == 0 {
            a0 = b0;
            ah = bh;
            z1 = z3;
            nz1 = nz3;
            continue;
        }
        // find which zeros of pl(x) are zeros of s(x).
        for i in 1..=j {
            if y[i - 1] < 0.0 || y[i - 1] > 1.0 {
                continue;
            }
            zero.push(t[l - 1] + h1 * y[i - 1]);
        }
        a0 = b0;
        ah = bh;
        z1 = z3;
        nz1 = nz3;
    }
    // the zeros of s(x) are arranged in increasing order.
    zero.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Array::from(zero)
}
