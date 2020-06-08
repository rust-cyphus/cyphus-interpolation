use ndarray::prelude::*;

pub(super) fn fpcurf(
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
    mut nrdata: ArrayViewMut1<usize>,
    ier: &mut i32,
) {
    let mut h = Array1::<f64>::zeros(7);
    let con1 = 0.1e0;
    let con9 = 0.9e0;
    let con4 = 0.4e-01;

    let mut acc = 0.0;
    let mut cos;
    let mut fpart;
    let mut fpms;
    let mut fpold = 0.0;
    let mut fp0 = 0.0;
    let mut f1;
    let mut f2;
    let mut f3;
    let mut p;
    let mut pinv;
    let mut piv;
    let mut p1;
    let mut p2;
    let mut p3;
    let mut rn;
    let mut sin;
    let mut store;
    let mut term;
    let mut wi;
    let mut xi;
    let mut yi;
    let mut i;
    let mut ich1;
    let mut ich3;
    let mut it: usize;
    let mut iter: usize;
    let mut i1;
    let mut i2;
    let mut i3;
    let mut j;
    let mut k3;
    let mut l;
    let mut l0;
    let mut mk1;
    let mut neww;
    let mut nk1;
    let mut nmax = 0;
    let mut nmin: usize;
    let mut nplus = 0;
    let mut npl1;
    let mut nrint;
    let mut n8;
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
    if iopt < 0 {
        // goto _60;
    }
    //  calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol * s;
    // determine nmax, the number of knots for spline interpolation.
    nmax = m + k1;
    if s > 0.0 {
        // goto _45;
    }
    // if s=0, s(x) is an interpolating spline.
    // test whether the required storage space exceeds the available 1.0.
    n = nmax;
    if nmax > nest {
        *ier = 1;
        return;
    }
    // find the position of the interior knots in case of interpolation.
    //_10:
    mk1 = m - k1;
    if mk1 == 0 {
        // goto _60;
    }
    k3 = k / 2;
    i = k2;
    j = k3 + 2;
    if k3 * 2 == k {
        // goto _30;
    }
    for l in 1..(mk1 + 1) {
        t[i - 1] = x[j - 1];
        i = i + 1;
        j = j + 1;
    }
    // goto _60;
    // _30:
    for l in 1..(mk1 + 1) {
        t[i - 1] = (x[j - 1] + x[j - 1 - 1]) * 0.5;
        i = i + 1;
        j = j + 1;
    }
    // goto _60;
    //   if s>0 our initial choice of knots depends on the value of iopt.
    //   if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
    //   polynomial of degree k which is a spline without interior knots.
    //   if iopt=1 and fp0>s we start computing the least squares spline
    //   according to the set of knots found at the last call of the routine.
    // _45:
    if iopt == 0 {
        // goto _50;
    }
    if n == nmin {
        // goto _50;
    }
    fp0 = fpint[n - 1];
    fpold = fpint[n - 1 - 1];
    nplus = nrdata[n - 1];
    if fp0 > s {
        // goto _60;
    }
    // _50:
    n = nmin;
    fpold = 0.0e0;
    nplus = 0;
    nrdata[0] = m - 2;
    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    // _60:
    for iter in 1..(m + 1) {
        if n == nmin {
            *ier = -2;
        }
        // find nrint, tne number of knot intervals.
        nrint = n - nmin + 1;
        // find the position of the additional knots which are needed for
        // the b-spline representation of s(x).
        nk1 = n - k1;
        i = n;
        for j in 1..(k1 + 1) {
            t[j - 1] = xb;
            t[i - 1] = xe;
            i = i - 1;
        }
        // compute the b-spline coefficients of the least-squares spline
        // sinf(x). the observation matrix a is built up row by row and
        // reduced to upper triangular form by givens transformations.
        // at the same time fp=f(p=inf) is computed.
        fp = 0.0;
        // initialize the observation matrix a.
        for i in 1..(nk1 + 1) {
            z[i - 1] = 0.0;
            for j in 1..(k1 + 1) {
                a[[i - 1, j - 1]] = 0.0;
            }
        }
        l = k1;
        for it in 1..(m + 1) {
            // fetch the current data point x(it),y(it).
            xi = x[it - 1];
            wi = w[it - 1];
            yi = y[it - 1] * wi;
            // search for knot interval t(l) <= xi < t(l+1).
            loop {
                if xi < t[l + 1 - 1] || l == nk1 {
                    break;
                }
                l = l + 1;
            }
            //  evaluate the (k+1) non-zero b-splines at xi and store them in q.
            // _90:
            super::fpbspl::fpbspl(t.view(), k, xi, l, h.view_mut());
            for i in 1..(k1 + 1) {
                q[[it - 1, i - 1]] = h[i - 1];
                h[i - 1] = h[i - 1] * wi;
            }
            // rotate the new row of the observation matrix into triangle.
            j = l - k1;
            for i in 1..(k1 + 1) {
                j = j + 1;
                piv = h[i - 1];
                if piv == 0.0 {
                    continue;
                }
                //  calculate the parameters of the givens transformation.
                super::fpgivs::fpgivs(piv, &mut a[[j - 1, 0]], cos, sin);
                // transformations to right hand side.
                super::fprota::fprota(*cos, *sin, &mut yi, &mut z[j - 1]);
                if i == k1 {
                    break;
                }
                i2 = 1;
                i3 = i + 1;
                for i1 in i3..(k1 + 1) {
                    i2 = i2 + 1;
                    // transformations to left hand side.
                    super::fprota::fprota(*cos, *sin, &mut h[i1 - 1], &mut a[[j - 1, i2 - 1]]);
                }
            }
            // add contribution of this row to the sum of squares of residual
            // right hand sides.
            // _120:
            fp = fp + yi * yi;
        }
        if *ier == -2 {
            fp0 = fp;
        }
        fpint[n - 1] = fp0;
        fpint[n - 1 - 1] = fpold;
        nrdata[n - 1] = nplus;
        // backward substitution to obtain the b-spline coefficients.
        super::fpback::fpback(a.view(), z.view(), nk1, k1, c);
        //  test whether the approximation sinf(x) is an acceptable solution.
        if iopt < 0 {
            return;
        }
        fpms = fp - s;
        if fpms.abs() < acc {
            return;
        }
        //  if f(p=inf) < s accept the choice of knots.
        if fpms < 0.0 {
            //goto _250;
        }
        //  if n = nmax, sinf(x) is an interpolating spline.
        if n == nmax {
            *ier = -1;
            return;
        }
        // increase the number of knots.
        // if n=nest we cannot increase the number of knots because of
        // the storage capacity limitation.
        if n == nest {
            *ier = 1;
            return;
        }
        //  determine the number of knots nplus we are going to add.
        if *ier == 0 {
            npl1 = nplus * 2;
            rn = nplus;
            if fpold - fp > acc {
                npl1 = (rn as f64 * fpms / (fpold - fp)) as usize;
            }
            nplus = (nplus * 2).max(npl1.max(nplus / 2).max(1));
        } else {
            nplus = 1;
            *ier = 0;
        }
        fpold = fp;
        // compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
        // t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.0;
        i = 1;
        l = k2;
        neww = 0;
        for it in 1..(m + 1) {
            if !(x[it - 1] < t[l - 1] || l > nk1) {
                neww = 1;
                l = l + 1;
            }
            term = 0.0;
            l0 = l - k2;
            for j in 1..(k1 + 1) {
                l0 = l0 + 1;
                term = term + c[l0 - 1] * q[[it - 1, j - 1]];
            }
            term = (w[it - 1] * (term - y[it - 1])).powi(2);
            fpart = fpart + term;
            if neww == 0 {
                continue;
            }
            store = term * 0.5;
            fpint[i - 1] = fpart - store;
            i = i + 1;
            fpart = store;
            neww = 0;
        }
        fpint[nrint - 1] = fpart;
        for l in 1..(nplus + 1) {
            // add a new knot.
            super::fpknot::fpknot(x.view(), t, &mut n, fpint, nrdata, &mut nrint, 1);
            // if n=nmax we locate the knots as for interpolation.
            if n == nmax {
                // goto _10;
            }
            // test whether we cannot further increase the number of knots.
            if n == nest {
                break;
            }
        }
        // restart the computations with the new set of knots.
        // _200:;
    }

    //   test whether the least-squares kth degree polynomial is a solution
    //   of our approximation problem.
    // _250:
    if *ier == -2 {
        return;
    }

    //
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
    //
    // evaluate the discontinuity jump of the kth derivative of the
    // b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
    super::fpdisc::fpdisc(t.view(), n, k2, b);
    // initial value for p.
    p1 = 0.0;
    f1 = fp0 - s;
    p3 = -1.0;
    f3 = fpms;
    p = 0.0;
    for i in 1..(nk1 + 1) {
        p = p + a[[i - 1, 0]];
    }
    rn = nk1;
    p = rn as f64 / p;
    ich1 = 0;
    ich3 = 0;
    n8 = n - nmin;
    //  iteration process to find the root of f(p) = s.
    for iter in 1..(maxit + 1) {
        // the rows of matrix b with weight 1/p are rotated into the
        // triangularised observation matrix a which is stored in g.
        pinv = 1.0 / p;
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
            yi = 0.0;
            for j in it..(nk1 + 1) {
                piv = h[0];
                // calculate the parameters of the givens transformation.
                super::fpgivs::fpgivs(piv, &mut g[[j - 1, 0]], cos, sin);
                // transformations to right hand side.
                super::fprota::fprota(*cos, *sin, &mut yi, &mut c[j - 1]);
                if j == nk1 {
                    // goto _300;
                }
                i2 = k1;
                if j > n8 {
                    i2 = nk1 - j;
                }
                for i in 1..(i2 + 1) {
                    //  transformations to left hand side.
                    i1 = i + 1;
                    super::fprota::fprota(*cos, *sin, &mut h[i1 - 1], &mut g[[j - 1, i1 - 1]]);
                    h[i - 1] = h[i1 - 1];
                }
                h[i2 + 1 - 1] = 0.0;
            }
            //_300:;
        }
        // backward substitution to obtain the b-spline coefficients.
        super::fpback::fpback(g.view(), c.view(), nk1, k2, c);
        // computation of f(p).
        fp = 0.0;
        l = k2;
        for it in 1..(m + 1) {
            if x[it - 1] < t[l - 1] || l > nk1 {
                //goto _310;
            }
            l = l + 1;
            // _310:
            l0 = l - k2;
            term = 0.0;
            for j in 1..(k1 + 1) {
                l0 = l0 + 1;
                term = term + c[l0 - 1] * q[[it - 1, j - 1]];
            }
            fp = fp + (w[it - 1] * (term - y[it - 1])).powi(2);
        }
        //  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp - s;
        if fpms.abs() < acc {
            return;
        }
        //  test whether the maximal number of iterations is reached.
        if iter == maxit {
            *ier = 3;
            return;
        }
        //  carry out 1.0 more step of the iteration process.
        p2 = p;
        f2 = fpms;
        if ich3 == 0 {
            if (f2 - f3) <= acc {
                //  our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if p <= p1 {
                    p = p1 * con9 + p2 * con1;
                }
                continue;
            }
            if f2 < 0.0 {
                ich3 = 1;
            }
        }
        if ich1 != 0 {
            continue;
        }
        if (f1 - f2) <= acc {
            //  our initial choice of p is too small
            p1 = p2;
            f1 = f2;
            p = p / con4;
            if p3 < 0.0 {
                continue;
            }
            if p >= p3 {
                p = p2 * con1 + p3 * con9;
            }
            continue;
        }
        if f2 > 0.0 {
            ich1 = 1;
        }
        //  test whether the iteration process proceeds as theoretically
        //  expected.
        // _350:
        if f2 >= f1 || f2 <= f3 {
            *ier = 2;
            return;
        }
        //  find the new value for p.
        p = super::fprati::fprati(&mut p1, &mut f1, p2, f2, &mut p3, &mut f3);
    }
    //   error codes and messages.
    // _400:
    // ier = 3;
    // return;
    // _410:
    // ier = 2;
    // return;
    // _420:
    // ier = 1;
    // return;
    // _430:
    // ier = -1;
}
