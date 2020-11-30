use ndarray::prelude::*;

/// Evaluates the (k+1) non-zero b-splines of degree k at t[l] <=x<=t[l+1]
/// using the stable recurrence relation of de Boor and Cox.
pub(super) fn fpbspl(t: ArrayView1<f64>, k: usize, x: f64, l: usize, mut h: ArrayViewMut1<f64>) {
    let mut hh = Array1::<f64>::zeros(19);
    h[0] = 1.0;
    for j in 1..(k + 1) {
        for i in 1..(j + 1) {
            hh[i - 1] = h[i - 1];
        }
        h[0] = 0.0;
        for i in 1..(j + 1) {
            let li = l + i;
            let lj = li - j;
            if t[li - 1] == t[lj - 1] {
                h[i + 1 - 1] = 0.0;
                continue;
            }
            let f = hh[i - 1] / (t[li - 1] - t[lj - 1]);
            h[i - 1] += f * (t[li - 1] - x);
            h[i + 1 - 1] = f * (x - t[lj - 1]);
        }
    }
}
