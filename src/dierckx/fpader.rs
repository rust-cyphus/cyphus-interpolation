use ndarray::prelude::*;

/// Calculates the derivatives of a spline of order k1 at the point
/// t[l] <= x < t[l+1], using the stable recurrence scheme of de Boor
pub(super) fn fpader(
    t: ArrayViewMut1<f64>,
    n: usize,
    c: ArrayViewMut1<f64>,
    k1: usize,
    x: f64,
    l: usize,
    d: ArrayViewMut1<f64>,
) {
    let mut h = Array1::<f64>::zeros(20);

    let lk = l - k1;
    for i in 1..(k1 + 1) {
        let ik = i + lk;
        h[i - 1] = c[ik - 1];
    }
    let mut kj = k1;
    let mut fac = 1.0;
    for j in 1..(k1 + 1) {
        let ki = kj;
        let j1 = j + 1;
        if j != 1 {
            let mut i = k1;
            for jj in j..(k1 + 1) {
                let li = i + lk;
                let lj = li + kj;
                h[i - 1] = (h[i - 1] - h[i - 1 - 1]) / (t[lj - 1] - t[li - 1]);
                i = i - 1;
            }
        }
        for i in j..(k1 + 1) {
            d[i - 1] = h[i - 1];
        }
        if j == k1 {
            for jj in j1..(k1 + 1) {
                let ki = ki - 1;
                let mut i = k1;
                for j2 in jj..(k1 + 1) {
                    let li = i + lk;
                    let lj = li + ki;
                    d[i - 1] = ((x - t[li - 1]) * d[i - 1] + (t[lj - 1] - x) * d[i - 1 - 1])
                        / (t[lj - 1] - t[li - 1]);
                    i = i - 1;
                }
            }
        }
        d[j - 1] = d[k1 - 1] * fac;
        let ak = k1 - j;
        fac = fac * ak;
        kj = kj - 1;
    }
}
