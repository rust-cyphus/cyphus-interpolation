use ndarray::prelude::*;

pub(super) fn splint(
    t: ArrayView1<f64>,
    n: usize,
    c: ArrayView1<f64>,
    k: usize,
    a: f64,
    b: f64,
) -> f64 {
    let nk1 = n - k - 1;
    let mut wrk = Array1::<f64>::zeros(nk1);
    // calculate the integrals wrk(i) of the normalized b-splines
    // ni,k+1(x), i=1,2,...nk1.
    super::fpintb::fpintb(t, n, wrk.view_mut(), nk1, a, b);
    // calculate the integral of s(x).
    c.iter()
        .take(nk1)
        .zip(wrk.iter().take(nk1))
        .fold(0.0, |acc, (z1, z2)| acc + (*z1) * (*z2))
}
