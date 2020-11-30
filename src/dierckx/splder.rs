use ndarray::prelude::*;

pub(super) fn splder(
    t: ArrayView1<f64>,
    n: usize,
    c: ArrayView1<f64>,
    k: usize,
    nu: usize,
    x: ArrayView1<f64>,
    mut y: ArrayViewMut1<f64>,
    m: usize,
    e: usize,
    mut wrk: ArrayViewMut1<f64>,
    ier: &mut usize,
) {
    let mut h = Array1::<f64>::zeros(6);
    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    // TODO: Removed nu < 0 check because I made nu a usize. Make sure this
    // is okay to do.
    if nu > k || m < 1 {
        return;
    }
    *ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    let k1 = k + 1;
    let k3 = k1 + 1;
    let nk1 = n - k1;
    let tb = t[k1 - 1];
    let te = t[nk1];
    let mut l1;
    // the derivative of order nu of a spline of degree k is a spline of
    // degree k-nu,the b-spline coefficients wrk(i) of which can be found
    // using the recurrence scheme of de boor.
    let mut l = 1;
    let mut kk = k;
    for i in 1..(nk1 + 1) {
        wrk[i - 1] = c[i - 1];
    }
    if nu != 0 {
        let mut nk2 = nk1;
        for _j in 1..(nu + 1) {
            let ak = kk;
            nk2 -= 1;
            l1 = l;
            for i in 1..(nk2 + 1) {
                l1 += 1;
                let l2 = l1 + kk;
                let fac = t[l2 - 1] - t[l1 - 1];
                if fac <= 0.0 {
                    continue;
                }
                wrk[i - 1] = ak as f64 * (wrk[i] - wrk[i - 1]) / fac;
            }
            l += 1;
            kk -= 1;
        }
        if kk == 0 {
            // if nu=k the derivative is a piecewise constant function
            let mut j = 1;
            for i in 1..(m + 1) {
                let arg = x[i - 1];
                // check if arg is in the support
                if arg < tb || arg > te {
                    if e == 1 {
                        y[i - 1] = 0.0;
                        continue;
                    } else if e == 2 {
                        *ier = 1;
                        return;
                    }
                }
                loop {
                    //  search for knot interval t(l) <= arg < t(l+1)
                    if arg >= t[l - 1] || l + 1 == k3 {
                        break;
                    }
                    // l1 = l;
                    l -= 1;
                    j -= 1;
                }

                loop {
                    if arg < t[l] || l == nk1 {
                        break;
                    }
                    l += 1;
                    j += 1;
                }
                y[i - 1] = wrk[j - 1];
            }
            return;
        }
    }
    l = k1;
    l1 = l + 1;
    let k2 = k1 - nu;
    //  main loop for the different points.
    for i in 1..(m + 1) {
        // fetch a new x-value arg.
        let arg = x[i - 1];
        // check if arg is in the support
        if arg < tb || arg > te {
            if e == 0 {
                // goto _135;
            } else if e == 1 {
                y[i - 1] = 0.0;
                continue;
            } else if e == 2 {
                *ier = 1;
                return;
            }
        }
        // search for knot interval t(l) <= arg < t(l+1)
        loop {
            if arg >= t[l - 1] || l1 == k3 {
                break;
            }
            l1 = l;
            l -= 1;
        }

        loop {
            if arg < t[l1 - 1] || l == nk1 {
                break;
            }
            l = l1;
            l1 = l + 1;
        }

        // evaluate the non-zero b-splines of degree k-nu at arg.
        // _150:
        super::fpbspl::fpbspl(t, kk, arg, l, h.view_mut());
        // find the value of the derivative at x=arg.
        let mut sp = 0.0;
        let mut ll = l - k1;
        for j in 1..(k2 + 1) {
            ll += 1;
            sp += wrk[ll - 1] * h[j - 1];
        }
        y[i - 1] = sp;
    }
}
