use ndarray::prelude::*;

#[allow(dead_code)]
pub(super) fn splev(
    t: ArrayView1<f64>,
    n: &mut usize,
    c: ArrayView1<f64>,
    k: &mut usize,
    x: ArrayView1<f64>,
    mut y: ArrayViewMut1<f64>,
    m: usize,
    e: usize,
    ier: &mut usize,
) {
    let mut h = Array1::<f64>::zeros(20);
    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    if m < 1 {
        return;
    }
    *ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    let k1 = *k + 1;
    let k2 = k1 + 1;

    let nk1 = *n - k1;
    let tb = t[k1 - 1];
    let te = t[nk1];
    let mut l = k1;
    let mut l1 = l + 1;

    // main loop for the different point
    for i in 1..(m + 1) {
        // fetch a new x-value arg.
        let mut arg = x[i - 1];
        if arg < tb || arg > te {
            if e == 1 {
                y[i - 1] = 0.0;
                continue;
            } else if e == 2 {
                *ier = 1;
                return;
            } else if e == 3 {
                if arg < tb {
                    arg = tb;
                } else {
                    arg = te;
                }
            }
        }
        // search for knot interval t(l) <= arg < t(l+1)
        loop {
            if arg >= t[l - 1] || l1 == k2 {
                break;
            }
            l1 = l;
            l = l - 1;
        }

        loop {
            if arg < t[l1 - 1] || l == nk1 {
                break;
            }
            l = l1;
            l1 = l + 1;
        }

        // evaluate the non-zero b-splines at arg.
        super::fpbspl::fpbspl(t, *k, arg, l, h.view_mut());
        //  find the value of s(x) at x=arg.
        let mut sp = 0.0;
        let mut ll = l - k1;
        for j in 1..(k1 + 1) {
            ll = ll + 1;
            sp = sp + c[ll - 1] * h[j - 1];
        }
        y[i - 1] = sp;
    }
}

#[allow(dead_code)]
pub(super) fn splev_point(
    t: ArrayView1<f64>,
    n: usize,
    c: ArrayView1<f64>,
    k: usize,
    x: f64,
    e: usize,
    ier: &mut usize,
) -> f64 {
    let mut y = 0.0;
    let mut h = Array1::<f64>::zeros(20);
    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    let k1 = k + 1;
    let k2 = k1 + 1;

    let nk1 = n - k1;
    let tb = t[k1 - 1];
    let te = t[nk1 + 1 - 1];
    let mut l = k1;
    let mut l1 = l + 1;

    // main loop for the different point
    // fetch a new x-value arg.
    let mut arg = x;
    if arg < tb || arg > te {
        if e == 1 {
            y = 0.0;
            return y;
        } else if e == 2 {
            *ier = 1;
            return y;
        } else if e == 3 {
            if arg < tb {
                arg = tb;
            } else {
                arg = te;
            }
        }
    }
    // search for knot interval t(l) <= arg < t(l+1)
    loop {
        if arg >= t[l - 1] || l1 == k2 {
            break;
        }
        l1 = l;
        l = l - 1;
    }

    loop {
        if arg < t[l1 - 1] || l == nk1 {
            break;
        }
        l = l1;
        l1 = l + 1;
    }

    // evaluate the non-zero b-splines at arg.
    super::fpbspl::fpbspl(t, k, arg, l, h.view_mut());
    //  find the value of s(x) at x=arg.
    let mut sp = 0.0;
    let mut ll = l - k1;
    for j in 1..(k1 + 1) {
        ll = ll + 1;
        sp = sp + c[ll - 1] * h[j - 1];
    }
    sp
}
