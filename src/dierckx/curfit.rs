use ndarray::prelude::*;

/// given the set of data points (x(i),y(i)) and the set of positive
/// numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
/// approximation of degree k on the interval xb <= x <= xe.
/// if iopt=-1 curfit calculates the weighted least-squares spline
/// according to a given set of knots.
/// if iopt>=0 the number of knots of the spline s(x) and the position
/// t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
/// ness of s(x) is then achieved by minimalizing the discontinuity
/// jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
/// n-k-1. the amount of smoothness is determined by the condition that
/// f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
/// negative constant, called the smoothing factor.
/// the fit s(x) is given in the b-spline representation (b-spline coef-
/// ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
/// subroutine splev.
#[allow(dead_code)]
pub(super) fn curfit(
    iopt: &mut i32,
    m: usize,
    x: ArrayView1<f64>,
    y: ArrayView1<f64>,
    w: ArrayView1<f64>,
    xb: f64,
    xe: f64,
    k: &mut usize,
    s: f64,
    nest: usize,
    n: &mut usize,
    mut t: ArrayViewMut1<f64>,
    mut c: ArrayViewMut1<f64>,
    fp: &mut f64,
    mut fpint: ArrayViewMut1<f64>,
    mut z: ArrayViewMut1<f64>,
    mut a: ArrayViewMut2<f64>,
    mut b: ArrayViewMut2<f64>,
    mut g: ArrayViewMut2<f64>,
    mut q: ArrayViewMut2<f64>,
    mut iwrk: ArrayViewMut1<usize>,
    ier: &mut i32,
) {
    // we set up the parameters tol and maxit
    let maxit = 20;
    let tol = 0.1e-02;
    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    *ier = 10;
    if *k <= 0 || *k > 5 {
        return;
    }
    let k1 = *k + 1;
    let k2 = k1 + 1;
    if *iopt < -1 || *iopt > 1 {
        return;
    }
    let nmin = 2 * k1;
    if m < k1 || nest < nmin {
        return;
    }
    //lwest = m * k1 + nest * (7 + 3 * k);
    //if (lwrk < lwest) return;
    if xb > x[0] || xe < x[m - 1] {
        return;
    }
    for i in 2..(m + 1) {
        if x[i - 2] > x[i - 1] {
            return;
        }
    }
    if *iopt < 0 {
        if *n < nmin || *n > nest {
            return;
        }
        let mut j = *n;
        for i in 1..(k1 + 1) {
            t[i - 1] = xb;
            t[j - 1] = xe;
            j -= 1;
        }
        *ier = super::fpchec::fpchec(x.view(), m, t.view_mut(), *n, *k);
        if *ier != 0 {
            return;
        }
    } else {
        if s < 0.0 {
            return;
        }
        if s == 0.0 && nest < (m + k1) {
            return;
        }
    }
    // we partition the working space and determine the spline approximation.
    super::fpcurf::fpcurf(
        *iopt,
        x.view(),
        y.view(),
        w.view(),
        m,
        xb,
        xe,
        k,
        s,
        nest,
        tol,
        maxit,
        k1,
        k2,
        n,
        t.view_mut(),
        c.view_mut(),
        fp,
        fpint.view_mut(),
        z.view_mut(),
        a.view_mut(),
        b.view_mut(),
        g.view_mut(),
        q.view_mut(),
        iwrk.view_mut(),
        ier,
    );
}
