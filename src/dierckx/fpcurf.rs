use ndarray::prelude::*;

/// Core routine of "curfit". Computes the knots 't' and spline coefficients 'c'
/// given x + y data, weights w, a degree k for the interpolating polynomials
/// and a smoothing factor s.
pub(super) fn fpcurf(
    iopt: i32,
    x: ArrayView1<f64>,
    y: ArrayView1<f64>,
    w: ArrayView1<f64>,
    m: usize,
    xb: f64,
    xe: f64,
    k: &mut usize,
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
    ier: &mut i32,
) {
    //double acc = 0, con1, con4, con9, cos, half, fpart, fpms, fpold = 0, fp0 = 0,
    // f1, f2, f3, one, p, pinv, piv, p1, p2, p3, rn, sin, store, term, wi,
    //xi, yi;
    //int i, ich1, ich3, it, iter, i1, i2, i3, j, k3, l, l0, mk1, neww, nk1,
    //nmax = 0, nmin, nplus = 0, npl1, nrint, n8;
    let mut h = Array1::<f64>::zeros(7);

    let con1 = 0.1e0;
    let con9 = 0.9e0;
    let con4 = 0.4e-01;
    let mut fpms = 0.0;
    let mut nk1 = 0;
    let mut nrint;

    let mut fp0 = 0.0;
    let mut fpold = 0.0;
    let mut nplus = 0;

    let mut goto60;

    let mut acc = 0.0;
    let mut nmax = 0;

    // part 1: determination of the number of knots and their position     c
    //  given a set of knots we compute the least-squares spline sinf(x),   c
    //  and the corresponding sum of squared residuals fp=f(p=inf).         c
    //  if iopt=-1 sinf(x) is the requested approximation.                  c
    //  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
    //  if fp <=s we will continue with the current set of knots.         c
    //  if fp > s we will increase the number of knots and compute the    c
    //      corresponding least-squares spline until finally fp<=s.        c
    //      the initial choice of knots depends on the value of s and iopt.   c
    //  if s=0 we have spline interpolation; in that case the number of   c
    //      knots equals nmax = m+k+1.                                        c
    //  if s > 0 and                                                      c
    //      iopt=0 we first compute the least-squares polynomial of         c
    //      degree k; n = nmin = 2*k+2                                      c
    //      iopt=1 we start with the set of knots found at the last         c
    //      call of the routine, except for the case that s > fp0; then     c
    //      we compute directly the least-squares polynomial of degree k.   c
    //  determine nmin, the number of knots for polynomial approximation.
    let nmin = 2 * k1;
    if iopt >= 0 {
        //  calculation of acc, the absolute tolerance for the root of f(p)=s.
        acc = tol * s;
        // determine nmax, the number of knots for spline interpolation.
        nmax = m + k1;
        if s > 0.0 {
            // if s>0 our initial choice of knots depends on the value of iopt.
            // if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
            // polynomial of degree k which is a spline without interior knots.
            // if iopt=1 and fp0>s we start computing the least squares spline
            // according to the set of knots found at the last call of the routine.
            //_45:
            if iopt == 0 || *n == nmin {
                *n = nmin;
                fpold = 0.0e0;
                nplus = 0;
                nrdata[0] = m - 2;
            } else {
                fp0 = fpint[*n - 1];
                fpold = fpint[*n - 1 - 1];
                nplus = nrdata[*n - 1];
                if fp0 <= s {
                    *n = nmin;
                    fpold = 0.0e0;
                    nplus = 0;
                    nrdata[0] = m - 2;
                }
            }
            goto60 = true;
        } else {
            // if s=0, s(x) is an interpolating spline.
            // test whether the required storage space exceeds the available one.
            *n = nmax;
            if nmax > nest {
                *ier = 1;
                return;
            } else {
                goto60 = false;
            }
        }
    } else {
        goto60 = true;
    }

    'find_knot: loop {
        if !goto60 {
            // find the position of the interior knots in case of interpolation.
            let mk1 = m - k1;
            if mk1 != 0 {
                let k3 = *k / 2;
                let mut i = k2;
                let mut j = k3 + 2;
                if k3 * 2 == *k {
                    for _l in 1..(mk1 + 1) {
                        t[i - 1] = (x[j - 1] + x[j - 1 - 1]) * 0.5;
                        i += 1;
                        j += 1;
                    }
                } else {
                    for _l in 1..(mk1 + 1) {
                        t[i - 1] = x[j - 1];
                        i += 1;
                        j += 1;
                    }
                }
            }
        }
        // main loop for the different sets of knots. m is a save upper bound
        // for the number of trials.
        for _iter in 1..(m + 1) {
            if *n == nmin {
                *ier = -2;
            }
            // find nrint, tne number of knot intervals.
            nrint = *n - nmin + 1;
            // find the position of the additional knots which are needed for
            // the b-spline representation of s(x).
            nk1 = *n - k1;
            let mut i = *n;
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
                super::fpbspl::fpbspl(t.view(), *k, xi, l, h.view_mut());
                for i in 1..(k1 + 1) {
                    q[[it - 1, i - 1]] = h[i - 1];
                    h[i - 1] *= wi;
                }
                // rotate the new row of the observation matrix into triangle.
                let mut j = l - k1;
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
            if *ier == -2 {
                fp0 = *fp;
            }
            fpint[*n - 1] = fp0;
            fpint[*n - 1 - 1] = fpold;
            nrdata[*n - 1] = nplus;
            // backward substitution to obtain the b-spline coefficients.
            super::fpback::fpback(a.view(), z.view(), nk1, k1, c.view_mut());
            //  test whether the approximation sinf(x) is an acceptable solution.
            if iopt < 0 {
                return;
            }
            fpms = *fp - s;
            if fpms.abs() < acc {
                return;
            }
            //  if f(p=inf) < s accept the choice of knots.
            if fpms < 0.0 {
                break 'find_knot;
            }
            //  if n = nmax, sinf(x) is an interpolating spline.
            if *n == nmax {
                *ier = -1;
                return;
            }
            // increase the number of knots.
            // if n=nest we cannot increase the number of knots because of
            // the storage capacity limitation.
            if *n == nest {
                *ier = 1;
                return;
            }
            //  determine the number of knots nplus we are going to add.
            if *ier == 0 {
                let mut npl1 = nplus * 2;
                let rn = nplus;
                if fpold - *fp > acc {
                    npl1 = rn * (fpms / (fpold - *fp)) as usize;
                }
                nplus = (nplus * 2).min(npl1.max(nplus / 2).max(1));
            } else {
                nplus = 1;
                *ier = 0;
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
                    l += 1;
                }
                let mut term = 0.0;
                let mut l0 = l - k2;
                for j in 1..(k1 + 1) {
                    l0 += 1;
                    term += c[l0 - 1] * q[[it - 1, j - 1]];
                }
                term = (w[it - 1] * (term - y[it - 1])).powi(2);
                fpart += term;
                if neww == 0 {
                    continue;
                }
                let store = term * 0.5;
                fpint[i - 1] = fpart - store;
                i += 1;
                fpart = store;
                neww = 0;
            }
            fpint[nrint - 1] = fpart;
            for _l in 1..(nplus + 1) {
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
                    goto60 = false;
                    continue 'find_knot;
                }
                // test whether we cannot further increase the number of knots.
                if *n == nest {
                    break;
                }
            }
            // restart the computations with the new set of knots.
        }
        break 'find_knot;
    }
    // test whether the least-squares kth degree polynomial is a solution
    // of our approximation problem.
    if *ier == -2 {
        return;
    }

    // part 2: determination of the smoothing spline sp(x).
    //
    // we have determined the number of knots and their position.
    // we now compute the b-spline coefficients of the smoothing spline
    // sp(x). the observation matrix a is extended by the rows of matrix
    // b expressing that the kth derivative discontinuities of sp(x) at
    // the interior knots t(k+2),...t(n-k-1) must be zero. the corres-
    // ponding weights of these additional rows are set to 1/p.
    // iteratively we then have to determine the value of p such that
    // f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that
    // the least-squares kth degree polynomial corresponds to p=0, and
    // that the least-squares spline corresponds to p=infinity. the
    // iteration process which is proposed here, makes use of rational
    // interpolation. since f(p) is a convex and strictly decreasing
    // function of p, it can be approximated by a rational function
    // r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-
    // ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used
    // to calculate the new value of p such that r(p)=s. convergence is
    // guaranteed by taking f1>0 and f3<0.
    // evaluate the discontinuity jump of the kth derivative of the
    // b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
    super::fpdisc::fpdisc(t.view(), *n, k2, b.view_mut());
    // initial value for p.
    let mut p1 = 0.0;
    let mut f1 = fp0 - s;
    let mut p3 = 0.0;
    let mut f3 = fpms;
    let mut p = 0.0;
    for i in 1..(nk1 + 1) {
        p += a[[i - 1, 0]];
    }
    let rn = nk1;
    p = rn as f64 / p;
    let mut ich1 = 0;
    let mut ich3 = 0;
    let n8 = *n - nmin;
    //  iteration process to find the root of f(p) = s.
    for iter in 1..(maxit + 1) {
        // the rows of matrix b with weight 1/p are rotated into the
        // triangularised observation matrix a which is stored in g.
        let pinv = 1.0 / p;
        for i in 1..(nk1 + 1) {
            c[i - 1] = z[i - 1];
            g[[i - 1, k2 - 1]] = 0.0;
            for j in 1..(k1 + 1) {
                g[[i - 1, j - 1]] = a[[i - 1, j - 1]];
            }
        }
        for it in 1..(n8 + 1) {
            // the row of matrix b is rotated into triangle by givens transformation
            for i in 1..(k2 + 1) {
                h[i - 1] = b[[it - 1, i - 1]] * pinv;
            }
            let mut yi = 0.0;
            for j in it..(nk1 + 1) {
                let piv = h[0];
                // calculate the parameters of the givens transformation.
                let mut cos = 0.0;
                let mut sin = 0.0;
                super::fpgivs::fpgivs(piv, &mut g[[j - 1, 0]], &mut cos, &mut sin);
                // transformations to right hand side.
                super::fprota::fprota(cos, sin, &mut yi, &mut c[j - 1]);
                if j == nk1 {
                    break;
                }
                let mut i2 = k1;
                if j > n8 {
                    i2 = nk1 - j;
                }
                for i in 1..(i2 + 1) {
                    //  transformations to left hand side.
                    let i1 = i + 1;
                    super::fprota::fprota(cos, sin, &mut h[i1 - 1], &mut g[[j - 1, i1 - 1]]);
                    h[i - 1] = h[i1 - 1];
                }
                h[i2 + 1 - 1] = 0.0;
            }
        }
        // backward substitution to obtain the b-spline coefficients.
        let mut cc = Array1::<f64>::zeros(c.len());
        cc.assign(&c);
        super::fpback::fpback(g.view(), cc.view(), nk1, k2, c.view_mut());
        // computation of f(p).
        *fp = 0.0;
        let mut l = k2;
        for it in 1..(m + 1) {
            if !(x[it - 1] < t[l - 1] || l > nk1) {
                l += 1;
            }
            let mut l0 = l - k2;
            let mut term = 0.0;
            for j in 1..(k1 + 1) {
                l0 += 1;
                term += c[l0 - 1] * q[[it - 1, j - 1]];
            }
            *fp += (w[it - 1] * (term - y[it - 1])).powi(2);
        }
        //  test whether the approximation sp(x) is an acceptable solution.
        fpms = *fp - s;
        if fpms.abs() < acc {
            return;
        }
        //  test whether the maximal number of iterations is reached.
        if iter == maxit {
            *ier = 3;
            return;
        }
        //  carry out one more step of the iteration process.
        let p2 = p;
        let f2 = fpms;
        if ich3 != 0 {
        } else if (f2 - f3) > acc {
            if f2 < 0.0 {
                ich3 = 1;
            }
        } else {
            //  our initial choice of p is too large.
            p3 = p2;
            f3 = f2;
            p *= con4;
            if p <= p1 {
                p = p1 * con9 + p2 * con1;
            }
            continue;
        }
        //_340:
        if ich1 != 0 {}
        if (f1 - f2) > acc {
            if f2 > 0.0 {
                ich1 = 1;
            }
        } else {
            //  our initial choice of p is too small
            p1 = p2;
            f1 = f2;
            p /= con4;
            if p3 < 0.0 {
                continue;
            }
            if p >= p3 {
                p = p2 * con1 + p3 * con9;
            }
            continue;
        }
        //  test whether the iteration process proceeds as theoretically
        //  expected.
        if f2 >= f1 || f2 <= f3 {
            *ier = 2;
            return;
        }
        //  find the new value for p.
        p = super::fprati::fprati(&mut p1, &mut f1, p2, f2, &mut p3, &mut f3);
    }
    *ier = 3;
}
