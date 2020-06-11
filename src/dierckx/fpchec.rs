use ndarray::prelude::*;

/// fpchec verifies the number and the position of the knots
/// t(j),j=1,2,...,n of a spline of degree k, in relation to the number
/// and the position of the data points x(i),i=1,2,...,m. if all of the
/// following conditions are fulfilled, the error parameter ier is set
/// to zero. if one of the conditions is violated ier is set to ten.
///
/// 1) k+1 <= n-k-1 <= m
/// 2) t(1) <= t(2) <= ... <= t(k+1)
///    t(n-k) <= t(n-k+1) <= ... <= t(n)
/// 3) t(k+1) < t(k+2) < ... < t(n-k)
/// 4) t(k+1) <= x(i) <= t(n-k)
/// 5) the conditions specified by schoenberg and whitney must hold
///    for at least one subset of data points, i.e. there must be a
///    subset of data points y(j) such that
///    t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
#[allow(dead_code)]
pub(super) fn fpchec(
    x: ArrayView1<f64>,
    m: usize,
    t: ArrayViewMut1<f64>,
    n: usize,
    k: usize,
) -> i32 {
    let mut ier: i32 = 10;
    let k1 = k + 1;
    let k2 = k1 + 1;
    let nk1 = n - k1;
    let nk2 = nk1 + 1;

    //  check condition no 1:
    // k+1 <= n-k-1 <= m
    if nk1 < k1 || nk1 > m {
        return ier;
    }
    // check condition no 2:
    // t(1) <= t(2) <= ... <= t(k+1)
    // t(n-k) <= t(n-k+1) <= ... <= t(n)
    let mut j = n;
    for i in 1..(k + 1) {
        if t[i - 1] > t[i] {
            return ier;
        }
        if t[j - 1] < t[j - 2] {
            return ier;
        }
        j = j - 1;
    }
    // check condition no 3:
    // t(k+1) < t(k+2) < ... < t(n-k)
    for i in k2..(nk2 + 1) {
        if t[i - 1] <= t[i - 2] {
            return ier;
        }
    }
    // check condition no 4:
    // t(k+1) <= x(i) <= t(n-k)
    if x[0] < t[k1 - 1] || x[m - 1] > t[nk2 - 1] {
        return ier;
    }
    // check condition no 5:
    // the conditions specified by schoenberg and whitney must hold
    // for at least one subset of data points, i.e. there must be a
    // subset of data points y(j) such that
    // t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
    if x[0] >= t[k2 - 1] || x[m - 1] <= t[nk1 - 1] {
        return ier;
    }
    let mut i = 1;
    let mut l = k2;
    let nk3 = nk1 - 1;
    if nk3 >= 2 {
        for j in 2..(nk3 + 1) {
            let tj = t[j - 1];
            l += 1;
            let tl = t[l - 1];
            loop {
                i += 1;
                if i >= m {
                    return ier;
                }
                if x[i - 1] <= tj {
                    continue;
                }
                if x[i - 1] >= tl {
                    return ier;
                }
                break;
            }
        }
    }
    ier = 0;
    ier
}
