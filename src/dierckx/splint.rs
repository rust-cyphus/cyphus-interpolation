use ndarray::prelude::*;

#[allow(dead_code)]
pub(super) fn splint(
    t: ArrayView1<f64>,
    n: usize,
    c: ArrayView1<f64>,
    k: usize,
    a: f64,
    b: f64,
    mut wrk: ArrayViewMut1<f64>,
) -> f64 {
    let nk1 = n - k - 1;
    // calculate the integrals wrk(i) of the normalized b-splines
    // ni,k+1(x), i=1,2,...nk1.
    super::fpintb::fpintb(t, n, wrk.view_mut(), nk1, a, b);
    // calculate the integral of s(x).
    let mut result = 0.0;
    for i in 1..(nk1 + 1) {
        result = result + c[i - 1] * wrk[i - 1];
    }
    result
}
