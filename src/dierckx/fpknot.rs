use ndarray::prelude::*;

/// locates an additional knot for a spline of degree
/// k and adjusts the corresponding parameters,i.e.
/// @param t The position of the knots.
/// @param n The number of knots.
/// @param fpint The sum of squares of residual right hand sides
///             for each knot interval.
/// @param nrdata The number of data points inside each knot interval.
/// @param nrint The number of knot intervals.
/// @param istart Indicates that the smallest data point at which the new knot
///              may be added is x(istart+1)
///
#[allow(dead_code)]
pub(super) fn fpknot(
    x: ArrayView1<f64>,
    mut t: ArrayViewMut1<f64>,
    n: &mut usize,
    mut fpint: ArrayViewMut1<f64>,
    mut nrdata: ArrayViewMut1<usize>,
    nrint: &mut usize,
    istart: usize,
) {
    let k = (*n - *nrint - 1) / 2;
    // search for knot interval t(number+k) <= x <= t(number+k+1) where
    // fpint(number) is maximal on the condition that nrdata(number)
    // not equals zero.
    let mut fpmax = 0.0;
    let mut jbegin = istart;
    let mut number = 0;
    let mut maxpt = 0;
    let mut maxbeg = 0;
    for j in 1..(*nrint + 1) {
        let jpoint = nrdata[j - 1];
        if !(fpmax >= fpint[j - 1] || jpoint == 0) {
            fpmax = fpint[j - 1];
            number = j;
            maxpt = jpoint;
            maxbeg = jbegin;
        }
        jbegin = jbegin + jpoint + 1;
    }
    // let coincide the new knot t(number+k+1) with a data point x(nrx)
    // inside the old knot interval t(number+k) <= x <= t(number+k+1).
    let ihalf = maxpt / 2 + 1;
    let nrx = maxbeg + ihalf;
    let next = number + 1;
    if next <= *nrint {
        // adjust the different parameters.
        for j in next..(*nrint + 1) {
            let jj = next + *nrint - j;
            fpint[jj] = fpint[jj - 1];
            nrdata[jj] = nrdata[jj - 1];
            let jk = jj + k;
            t[jk] = t[jk - 1];
        }
    }
    nrdata[number - 1] = ihalf - 1;
    nrdata[next - 1] = maxpt - ihalf;
    let am = maxpt;
    let mut an = nrdata[number - 1];
    fpint[number - 1] = fpmax * an as f64 / am as f64;
    an = nrdata[next - 1];
    fpint[next - 1] = fpmax * an as f64 / am as f64;
    let jk = next + k;
    t[jk - 1] = x[nrx - 1];
    *n += 1;
    *nrint += 1;
}
