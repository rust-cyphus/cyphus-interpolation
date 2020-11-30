use ndarray::prelude::*;

/// calculates integrals of the normalized b-splines
/// nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
/// it makes use of the formulae of gaffney for the calculation of
/// indefinite integrals of b-splines.
///
/// @param t Length n, containing the position of the knots.
/// @param n integer value, giving the number of knots.
/// @param bint (out) array,length nk1, containing the integrals of the
/// b-splines.
/// @param nk1 integer value, giving the number of b-splines of degree k,
///            defined on the set of knots ,i.e. nk1 = n-k-1.
/// @param x lower bound
/// @param y upper bound
pub(super) fn fpintb(
    t: ArrayView1<f64>,
    n: usize,
    mut bint: ArrayViewMut1<f64>,
    nk1: usize,
    x: f64,
    y: f64,
) {
    let mut aint = Array1::<f64>::zeros(6);
    let mut h = Array1::<f64>::zeros(6);
    let mut h1 = Array1::<f64>::zeros(6);

    // initialization.
    let k1 = n - nk1;
    let ak = k1;
    let k = k1 - 1;
    bint.fill(0.0);
    // the integration limits are arranged in increasing order.
    let mut a = x;
    let mut b = y;
    let mut min = 0;
    if a >= b {
        if a == b {
            return;
        }
        a = y;
        b = x;
        min = 1;
    }

    if a < t[k1 - 1] {
        a = t[k1 - 1];
    }
    if b > t[nk1] {
        b = t[nk1];
    }
    if a > b {
        return;
    }

    // using the expression of gaffney for the indefinite integral of a
    // b-spline we find that
    // bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
    //     where for t(l) <= x < t(l+1)
    //     res(j,x) = 0, j=1,2,...,l-k-1
    //              = 1, j=l+1,l+2,...,nk1
    //              = aint(j+k-l+1), j=l-k,l-k+1,...,l
    //                = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
    //                  i=0,1,...,k
    let mut l = k1;
    let mut l0 = l + 1;
    let mut arg = a;
    let mut ia = 0;
    for it in 1..=2 {
        //  search for the knot interval t(l) <= arg < t(l+1).
        loop {
            if arg < t[l0 - 1] || l == nk1 {
                break;
            }
            l = l0;
            l0 += 1;
        }
        // calculation of aint(j), j=1,2,...,k+1.
        // initialization.
        for j in 0..=k1 {
            aint[j] = 0.0;
        }
        aint[0] = (arg - t[l - 1]) / (t[l] - t[l - 1]);
        h1[0] = 1.0;
        for j in 1..=k {
            // evaluation of the non-zero b-splines of degree j at arg,i.e.
            // h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
            h[0] = 0.0;
            for i in 1..=j {
                let li = l + i;
                let lj = li - j;
                let f = h1[i - 1] / (t[li - 1] - t[lj - 1]);
                h[i - 1] += f * (t[li - 1] - arg);
                h[i] = f * (arg - t[lj - 1]);
            }
            // updating of the integrals aint.
            let j1 = j + 1;
            for i in 1..=j1 {
                let li = l + i;
                let lj = li - j1;
                aint[i - 1] += h[i - 1] * (arg - t[lj - 1]) / (t[li - 1] - t[lj - 1]);
                h1[i - 1] = h[i - 1];
            }
        }
        if it == 1 {
            //  updating of the integrals bint
            let mut lk = l - k;
            ia = lk;
            for i in 1..=k1 {
                bint[lk - 1] = -aint[i - 1];
                lk += 1;
            }
            arg = b;
        }
    }
    //  updating of the integrals bint.
    let mut lk = l - k;
    let ib = lk - 1;
    for i in 1..=k1 {
        bint[lk - 1] += aint[i - 1];
        lk += 1;
    }
    if ib >= ia {
        for i in ia..=ib {
            bint[i - 1] += 1.0;
        }
    }
    // the scaling factors are taken into account.
    let f = 1.0 / ak as f64;
    for i in 1..=nk1 {
        let j = i + k1;
        bint[i - 1] *= (t[j - 1] - t[i - 1]) * f;
    }
    // the order of the integration limits is taken into account.
    if min == 0 {
        return;
    }
    for i in 0..nk1 {
        bint[i] = -bint[i];
    }
}
