use ndarray::prelude::*;

/// Find the real zero of a linear polynomial
fn real_zeros_linear(a: f64, b: f64, c: f64, d: f64, mut x: ArrayViewMut1<f64>) -> usize {
    x[0] = -d / c;
    apply_newton_iter(a, b, c, d, x, 1)
}

/// Find the real zeros of a quadratic polynomial.
fn real_zeros_quadratic(a: f64, b: f64, c: f64, d: f64, mut x: ArrayViewMut1<f64>) -> usize {
    let disc = c * c - 4.0 * b * d;
    if disc < 0.0 {
        return 0;
    }
    let u = disc.sqrt();
    let b1 = b + b;
    x[0] = (-c + u) / b1;
    x[1] = (-c - u) / b1;

    apply_newton_iter(a, b, c, d, x, 2)
}

/// Find the real zeros of a cubic polynomial.
fn real_zeros_cubic(a: f64, b: f64, c: f64, d: f64, mut x: ArrayViewMut1<f64>) -> usize {
    let pi3 = std::f64::consts::PI / 3.0;
    let b1 = b / a * (1.0 / 3.0);
    let c1 = c / a;
    let d1 = d / a;
    let q = c1 * (1.0 / 3.0) - b1 * b1;
    let r = b1 * b1 * b1 + (d1 - b1 * c1) * 0.5;
    let disc = q * q * q + r * r;
    if disc <= 0.0 {
        let mut u = q.abs().sqrt();
        if r < 0.0 {
            u = -u;
        }
        let p3 = (-disc).sqrt().atan2(r.abs() * (1.0 / 3.0));
        let u2 = u + u;
        x[0] = -u2 * p3.cos() - b.abs();
        x[1] = u2 * (pi3 - p3).cos() - b.abs();
        x[2] = u2 * (pi3 + p3).cos() - b.abs();
        return apply_newton_iter(a, b, c, d, x, 3);
    }
    let u = disc.sqrt();
    let u1 = -r + u;
    let u2 = -r - u;
    x[0] = u1.abs().powf(1.0 / 3.0).copysign(u1) + u2.abs().powf(1.0 / 3.0).copysign(u2) - b.abs();
    apply_newton_iter(a, b, c, d, x, 1)
}

/// Apply a few newton iterations to improve the accuracy of the roots.
fn apply_newton_iter(a: f64, b: f64, c: f64, d: f64, mut x: ArrayViewMut1<f64>, n: usize) -> usize {
    for i in 1..(n + 1) {
        let y = x[i - 1];
        let f = ((a * y + b) * y + c) * y + d;
        let df = (3.0 * a * y + 2.0 * b) * y + c;
        let mut step = 0.0;
        if f.abs() < df.abs() * 0.1 {
            step = f / df;
        }
        x[i - 1] = y - step;
    }
    n
}

/// Find the real zeros of a cubic polynomial p(x) = ax^3 + bx^2 + cx + d
pub(super) fn fpcuro(a: f64, b: f64, c: f64, d: f64, mut x: ArrayViewMut1<f64>) -> usize {
    // test whether p(x) is a third degree polynomial.
    if b.abs().max(c.abs()).max(d.abs()) >= a.abs() * 1e4 {
        // test whether p(x) is a second degree polynomial.
        if c.abs().max(d.abs()) >= b.abs() * 1e4 {
            // test whether p(x) is a first degree polynomial.
            if d.abs() >= c.abs() * 1e4 {
                // p(x) is a constant function.
                return 0;
            }
            //  p(x) is a first degree polynomial.
            return real_zeros_linear(a, b, c, d, x);
        }
        //  p(x) is a second degree polynomial.
        return real_zeros_quadratic(a, b, c, d, x);
    }
    real_zeros_cubic(a, b, c, d, x)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_linear() {
        let (a, b, c, d) = (0.0, 0.0, 1.0, 1.0);
        let mut x = array![0.0, 0.0, 0.0, 0.0];
        let n = fpcuro(a, b, c, d, x.view_mut());
        assert!((x[0] + 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_quadratic() {
        let (a, b, c, d) = (0.0, 1.0, 0.0, -2f64);
        let mut x = array![0.0, 0.0, 0.0, 0.0];
        let n = fpcuro(a, b, c, d, x.view_mut());
        assert!((x[0] - 2f64.sqrt()).abs() < 1e-5);
        assert!((x[1] + 2f64.sqrt()).abs() < 1e-5);
    }
}
