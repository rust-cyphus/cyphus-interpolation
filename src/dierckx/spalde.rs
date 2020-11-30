use ndarray::prelude::*;

pub(super) fn spalde(
    t: ArrayView1<f64>,
    n: usize,
    c: ArrayView1<f64>,
    k1: usize,
    x: f64,
    d: ArrayViewMut1<f64>,
    ier: &mut usize,
) {
    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    let nk1 = n - k1;
    if x < t[k1 - 1] || x > t[nk1 + 1 - 1] {
        return;
    }
    // search for knot interval t(l) <= x < t(l+1)
    let mut l = k1;
    loop {
        if x < t[l + 2 - 1] || l == nk1 {
            break;
        }
        l += 1;
    }
    if t[l - 1] >= t[l + 1 - 1] {
        return;
    }
    *ier = 0;
    //  calculate the derivatives.
    super::fpader::fpader(t, c, k1, x, l, d);
}
