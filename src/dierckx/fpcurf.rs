use ndarray::prelude::*;

/// Determination of the number of knots and their position
/// given a set of knots we compute the least-squares spline sinf(x),
/// and the corresponding sum of squared residuals fp=f(p=inf).
/// if iopt=-1 sinf(x) is the requested approximation.
/// if iopt=0 or iopt=1 we check whether we can accept the knots:
/// if fp <=s we will continue with the current set of knots.
/// if fp > s we will increase the number of knots and compute the
///     corresponding least-squares spline until finally fp<=s.
///     the initial choice of knots depends on the value of s and iopt.
/// if s=0 we have spline interpolation; in that case the number of
///     knots equals nmax = m+k+1.
/// if s > 0 and
///     iopt=0 we first compute the least-squares polynomial of
///     degree k; n = nmin = 2*k+2
///     iopt=1 we start with the set of knots found at the last
///     call of the routine, except for the case that s > fp0; then
///     we compute directly the least-squares polynomial of degree k.
/// determine nmin, the number of knots for polynomial approximation.
pub(super) fn construct_knots(
    iopt: usize,
    mut x: ArrayViewMut1<f64>,
    mut y: ArrayViewMut1<f64>,
    mut w: ArrayViewMut1<f64>,
    m: usize,
    xb: f64,
    xe: f64,
    k: usize,
    s: f64,
    nest: usize,
    tol: f64,
    maxit: usize,
    k1: usize,
    k2: usize,
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
    mut nrdata: ArrayViewMut1<usize>,
) -> i32 {
    let mut h = Array1::<f64>::zeros(7);
    let mut ier: i32 = 0;

    let nmin = 2 * k1;
    if iopt < 0 {
        // goto _60;
    }
    //  calculation of acc, the absolute tolerance for the root of f(p)=s.
    let acc = tol * s;
    // determine nmax, the number of knots for spline interpolation.
    let nmax = m + k1;
    if s > 0.0 {
        //goto _45;
    }
    // if s=0, s(x) is an interpolating spline.
    // test whether the required storage space exceeds the available one.
    *n = nmax;
    if nmax > nest {
        //goto _420;
    }
    // find the position of the interior knots in case of interpolation.
    //_10:
    let mk1 = m - k1;
    if mk1 == 0 {
        //goto _60;
    }
    let k3 = k / 2;
    let mut i = k2;
    let mut j = k3 + 2;
    if k3 * 2 == k {
        //goto _30;
    }
    for l in 1..(mk1 + 1) {
        t[i - 1] = x[j - 1];
        i += 1;
        j += 1;
    }
    //goto _60;
    //_30:
    for l in 1..(mk1 + 1) {
        t[i - 1] = (x[j - 1] + x[j - 1 - 1]) * 0.5;
        i += 1;
        j += 1;
    }
    //goto _60;
    // if s>0 our initial choice of knots depends on the value of iopt.
    // if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
    // polynomial of degree k which is a spline without interior knots.
    // if iopt=1 and fp0>s we start computing the least squares spline
    // according to the set of knots found at the last call of the routine.
    //_45:
    if iopt == 0 {
        //goto _50;
    }
    if *n == nmin {
        //goto _50;
    }
    let mut fp0 = fpint[*n - 1];
    let mut fpold = fpint[*n - 1 - 1];
    let mut nplus = nrdata[*n - 1];
    if fp0 > s {
        //goto _60;
    }
    //_50:
    *n = nmin;
    fpold = 0.0e0;
    nplus = 0;
    nrdata[0] = m - 2;
    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    //_60:
    for _ in 1..(m + 1) {
        if *n == nmin {
            ier = -2;
        }
        // find nrint, tne number of knot intervals.
        let mut nrint = *n - nmin + 1;
        // find the position of the additional knots which are needed for
        // the b-spline representation of s(x).
        let nk1 = *n - k1;
        i = *n;
        for j in 1..(k1 + 1) {
            t[j - 1] = xb;
            t[i - 1] = xe;
            i -= 1;
        }
        // compute the b-spline coefficients of the least-squares spline
        // sinf(x). the observation matrix a is built up row by row and
        // reduced to upper triangular form by givens transformations.
        // at the same time fp=f(p=inf) is computed.
        *fp = 0.0;
        // initialize the observation matrix a.
        for i in 1..(nk1 + 1) {
            z[i - 1] = 0.0;
            for j in 1..(k1 + 1) {
                a[[i - 1, j - 1]] = 0.0;
            }
        }
        let mut l = k1;
        for it in 1..(m + 1) {
            // fetch the current data point x(it),y(it).
            let xi = x[it - 1];
            let wi = w[it - 1];
            let mut yi = y[it - 1] * wi;
            // search for knot interval t(l) <= xi < t(l+1).
            loop {
                if xi < t[l + 1 - 1] || l == nk1 {
                    break;
                }
                l += 1;
            }
            //  evaluate the (k+1) non-zero b-splines at xi and store them in q.
            super::fpbspl::fpbspl(t.view(), k, xi, l, h.view_mut());
            for i in 1..(k1 + 1) {
                q[[it - 1, i - 1]] = h[i - 1];
                h[i - 1] *= wi;
            }
            // rotate the new row of the observation matrix into triangle.
            j = l - k1;
            for i in 1..(k1 + 1) {
                j += 1;
                let piv = h[i - 1];
                if piv == 0.0 {
                    continue;
                }
                //  calculate the parameters of the givens transformation.
                let mut cos = 0.0;
                let mut sin = 0.0;
                super::fpgivs::fpgivs(piv, &mut a[[j - 1, 0]], &mut cos, &mut sin);
                // transformations to right hand side.
                super::fprota::fprota(cos, sin, &mut yi, &mut z[j - 1]);
                if i == k1 {
                    break;
                }
                let mut i2 = 1;
                let i3 = i + 1;
                for i1 in i3..(k1 + 1) {
                    i2 += 1;
                    // transformations to left hand side.
                    super::fprota::fprota(cos, sin, &mut h[i1 - 1], &mut a[[j - 1, i2 - 1]]);
                }
            }
            // add contribution of this row to the sum of squares of residual
            // right hand sides.
            *fp += yi * yi;
        }
        if ier == -2 {
            fp0 = *fp;
        }
        fpint[*n - 1] = fp0;
        fpint[*n - 1 - 1] = fpold;
        nrdata[*n - 1] = nplus;
        // backward substitution to obtain the b-spline coefficients.
        super::fpback::fpback(a.view(), z.view(), nk1, k1, c.view_mut());
        //  test whether the approximation sinf(x) is an acceptable solution.
        if iopt < 0 {
            return ier;
        }
        let fpms = *fp - s;
        if fpms.abs() < acc {
            return ier;
        }
        //  if f(p=inf) < s accept the choice of knots.
        if fpms < 0.0 {
            //goto _250;
        }
        //  if n = nmax, sinf(x) is an interpolating spline.
        if *n == nmax {
            //goto _430;
        }
        // increase the number of knots.
        // if n=nest we cannot increase the number of knots because of
        // the storage capacity limitation.
        if *n == nest {
            //goto _420;
        }
        //  determine the number of knots nplus we are going to add.
        if ier == 0 {
            //goto _140;
            let mut npl1 = nplus * 2;
            let rn = nplus;
            if fpold - *fp > acc {
                npl1 = rn * (fpms / (fpold - *fp)) as usize;
            }
            nplus = (nplus * 2).min(npl1.max(nplus / 2).max(1));
        } else {
            nplus = 1;
            ier = 0;
        }
        fpold = *fp;
        // compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
        // t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        let mut fpart = 0.0;
        i = 1;
        l = k2;
        let mut neww = 0;
        for it in 1..(m + 1) {
            if !(x[it - 1] < t[l - 1] || l > nk1) {
                neww = 1;
                l = l + 1;
            }
            let mut term = 0.0;
            let mut l0 = l - k2;
            for j in 1..(k1 + 1) {
                l0 = l0 + 1;
                term += c[l0 - 1] * q[[it - 1, j - 1]];
            }
            term = (w[it - 1] * (term - y[it - 1])).powi(2);
            fpart = fpart + term;
            if neww == 0 {
                continue;
            }
            let store = term * 0.5;
            fpint[i - 1] = fpart - store;
            i = i + 1;
            fpart = store;
            neww = 0;
        }
        fpint[nrint - 1] = fpart;
        for _ in 1..(nplus + 1) {
            // add a new knot.
            super::fpknot::fpknot(
                x.view(),
                t.view_mut(),
                n,
                fpint.view_mut(),
                nrdata.view_mut(),
                &mut nrint,
                1,
            );
            // if n=nmax we locate the knots as for interpolation.
            if *n == nmax {
                //goto _10;
            }
            // test whether we cannot further increase the number of knots.
            if *n == nest {
                break;
            }
        }
        // restart the computations with the new set of knots.
    }
    // test whether the least-squares kth degree polynomial is a solution
    // of our approximation problem.
    //_250:
    if ier == -2 {
        return ier;
    }
    ier
}
