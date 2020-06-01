use ndarray::prelude::*;

pub(super) fn fpdisc(t: ArrayView1<f64>, n: usize, k2: usize, mut b: ArrayViewMut2<f64>) {
    let mut h = Array1::<f64>::zeros(12);
    let k1 = k2 - 1;
    let k = k1 - 1;
    let nk1 = n - k1;
    let nrint = nk1 - k;
    let an = nrint;
    let fac = an as f64 / (t[nk1] - t[k1 - 1]);
    for l in k2..(nk1 + 1) {
        let lmk = l - k1;
        for j in 1..(k1 + 1) {
            let ik = j + k1;
            let lj = l + j;
            let lk = lj - k2;
            h[j - 1] = t[l - 1] - t[lk - 1];
            h[ik - 1] = t[l - 1] - t[lj - 1];
        }
        let mut lp = lmk;
        for j in 1..(k2 + 1) {
            let mut jk = j;
            let mut prod = h[j - 1];
            for _i in 1..(k + 1) {
                jk = jk + 1;
                prod = prod * h[jk - 1] * fac;
            }
            let lk = lp + k1;
            b[[lmk - 1, j - 1]] = (t[lk - 1] - t[lp - 1]) / prod;
            lp = lp + 1;
        }
    }
}
