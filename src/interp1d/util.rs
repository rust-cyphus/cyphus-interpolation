use super::acc::InterpAccel;

#[allow(dead_code)]
pub(super) fn bsearch(xarr: &[f64], x: f64, idx_low: usize, idx_high: usize) -> usize {
    let mut ilow = idx_low;
    let mut ihigh = idx_high;

    while ihigh > ilow + 1 {
        let i = (ihigh + ilow) / 2;
        if xarr[i] > x {
            ihigh = i;
        } else {
            ilow = i;
        }
    }
    ilow
}

#[allow(dead_code)]
pub(super) fn accel_find(xarr: &[f64], x: f64, acc: &mut InterpAccel) -> usize {
    let xidx = acc.cache;

    if x < xarr[xidx] {
        acc.miss_count += 1;
        acc.cache = bsearch(xarr, x, 0, xidx);
    } else if x >= xarr[xidx + 1] {
        acc.miss_count += 1;
        acc.cache = bsearch(xarr, x, xidx, xarr.len() - 1);
    } else {
        acc.hit_count += 1;
    }
    acc.cache
}

#[allow(dead_code)]
pub(super) fn integ_eval(ai: f64, bi: f64, ci: f64, di: f64, x: f64, a: f64, b: f64) -> f64 {
    let r1 = a - x;
    let r2 = b - x;
    let r12 = r1 + r2;
    let bterm = 0.5 * bi * r12;
    let cterm = (1.0 / 3.0) * ci * (r1 * r1 + r2 * r2 + r1 * r2);
    let dterm = 0.25 * di * r12 * (r1 * r1 + r2 * r2);

    (b - a) * (ai + bterm + cterm + dterm)
}

#[allow(dead_code)]
pub(super) fn solve_tridiag(
    diag: &[f64],
    offdiag: &[f64],
    b: &[f64],
    x: &mut [f64],
    n: usize,
) -> bool {
    let mut gamma = vec![0.0; n];
    let mut alpha = vec![0.0; n];
    let mut c = vec![0.0; n];
    let mut z = vec![0.0; n];

    alpha[0] = diag[0];
    gamma[0] = offdiag[0] / alpha[0];

    if alpha[0] == 0.0 {
        return false;
    }

    for i in 1..(n - 1) {
        alpha[i] = diag[i] - offdiag[i - 1] * gamma[i - 1];
        gamma[i] = offdiag[i] / alpha[i];
        if alpha[i] == 0.0 {
            return false;
        }
    }

    if n > 1 {
        alpha[n - 1] = diag[n - 1] - offdiag[n - 2] * gamma[n - 2];
    }

    // update RHS
    z[0] = b[0];
    for i in 1..n {
        z[i] = b[i] - gamma[i - 1] * z[i - 1];
    }
    for i in 0..n {
        c[i] = z[i] / alpha[i];
    }

    // backsubstitution
    x[n - 1] = c[n - 1];
    if n >= 2 {
        for i in (0..(n - 1)).rev() {
            x[i] = c[i] - gamma[i] * x[i + 1];
        }
    }
    true
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bsearch() {
        let xarr = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        assert_eq!(bsearch(&xarr, 1.5, 0, 4), 1);
        assert_eq!(bsearch(&xarr, 4.0, 0, 4), 3);
        assert_eq!(bsearch(&xarr, 0.0, 0, 4), 0);
        assert_eq!(bsearch(&xarr, 2.0, 0, 4), 2);
        assert_eq!(bsearch(&xarr, 10.0, 0, 4), 3);
        assert_eq!(bsearch(&xarr, -10.0, 0, 4), 0);
    }
    #[test]
    fn test_accel() {
        let xarr = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let mut k1 = 0;
        let mut k2 = 0;
        let mut t = false;

        let r = vec![
            -0.2, 0.0, 0.1, 0.7, 1.0, 1.3, 1.9, 2.0, 2.2, 2.7, 3.0, 3.1, 3.6, 4.0, 4.1, 4.9,
        ];
        let mut acc = InterpAccel::new();

        while k1 < 16 && k2 < 16 {
            let x = r[if t { k1 } else { k2 }];
            t = !t;
            if !t {
                k1 = (k1 + 1) % 16;
                if k1 == 0 {
                    k2 += 1;
                }
            }

            let i = accel_find(&xarr, x, &mut acc);
            let j = bsearch(&xarr, x, 0, 4);
            assert_eq!(i, j);
        }
    }
}
