use ndarray::prelude::*;

pub(super) fn curfit(
    iopt: usize,
    m: usize,
    mut x: ArrayViewMut1<f64>,
    mut y: ArrayViewMut1<f64>,
    mut w: ArrayViewMut1<f64>,
    xb: f64,
    xe: f64,
    k: usize,
    s: f64,
    nest: usize,
    n: usize,
    mut t: ArrayViewMut1<f64>,
    mut c: ArrayViewMut1<f64>,
    fp: f64,
    mut fpint: ArrayViewMut1<f64>,
    mut z: ArrayViewMut1<f64>,
    mut a: ArrayViewMut2<f64>,
    mut b: ArrayViewMut2<f64>,
    mut g: ArrayViewMut2<f64>,
    mut q: ArrayViewMut2<f64>,
    mut iwrk: ArrayViewMut1<usize>,
    ier: &mut usize,
) {
    // we set up the parameters tol and maxit
    maxit = 20;
    tol = 0.1e-02;
    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    ier = 10;
    if (k <= 0 || k > 5) return;
    k1 = k + 1;
    k2 = k1 + 1;
    if ((iopt < (-1)) || iopt > 1) return;
    nmin = 2 * k1;
    if (m < k1 || nest < nmin) return;
    //lwest = m * k1 + nest * (7 + 3 * k);
    //if (lwrk < lwest) return;
    if (xb > x[0] || xe < x[m - 1]) return;
    for (i = 2; i <= m; i++) {
        if (x[i - 2] > x[i - 1]) return;
    }
    if (iopt >= 0) goto _30;
    if (n < nmin || n > nest) return;
    j = n;
    for (i = 1; i <= k1; i++) {
        t[i - 1] = xb;
        t[j - 1] = xe;
        j = j - 1;
    }
    fpchec(x, m, t, n, k, ier);
    if (ier == 0) goto _40;
    return;
    _30:
    if (s < 0.0) return;
    if (s == 0.0 && nest < (m + k1)) return;
    // we partition the working space and determine the spline approximation.
    _40:

    fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2, n, t, c, fp,
           fpint, z, a, b, g, q, iwrk, ier);
}
