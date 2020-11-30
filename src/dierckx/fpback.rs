use ndarray::prelude::*;

/// Calculate the solution of the system of equations `a`*`c` = `z` where `a` is
/// an `n`x`n` upper triangular matrix with bandwidth `k`.
pub(super) fn fpback(
    a: ArrayView2<f64>,
    z: ArrayView1<f64>,
    n: usize,
    k: usize,
    mut c: ArrayViewMut1<f64>,
) {
    let k1 = k - 1;
    c[n - 1] = z[n - 1] / a[[n - 1, 0]];
    let mut i = n - 1;
    if i == 0 {
        return;
    }

    for j in 2..(n + 1) {
        let mut store = z[i - 1];
        let mut i1 = k1;
        if j <= k1 {
            i1 = j - 1;
        }
        let mut m = i;
        for l in 1..(i1 + 1) {
            m += 1;
            store -= c[m - 1] * a[[i - 1, l + 1 - 1]];
        }
        c[i - 1] = store / a[[i - 1, 0]];
        i -= 1;
    }
}
