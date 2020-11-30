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
        let u = (if r < 0.0 { -1.0 } else { 1.0 }) * q.abs().sqrt();
        let p3 = (-disc).sqrt().atan2(r.abs()) * (1.0 / 3.0);
        let u2 = u + u;
        x[0] = -u2 * p3.cos() - b1;
        x[1] = u2 * (pi3 - p3).cos() - b1;
        x[2] = u2 * (pi3 + p3).cos() - b1;
        return apply_newton_iter(a, b, c, d, x, 3);
    }
    let u = disc.sqrt();
    let u1 = -r + u;
    let u2 = -r - u;
    x[0] = u1.abs().powf(1.0 / 3.0).copysign(u1) + u2.abs().powf(1.0 / 3.0).copysign(u2) - b1;
    apply_newton_iter(a, b, c, d, x, 1)
}

/// Apply a few newton iterations to improve the accuracy of the roots.
fn apply_newton_iter(a: f64, b: f64, c: f64, d: f64, mut x: ArrayViewMut1<f64>, n: usize) -> usize {
    for i in 0..n {
        let y = x[i];
        let f = a.mul_add(y, b).mul_add(y, c).mul_add(y, d);
        let df = (3.0 * a).mul_add(y, 2.0 * b).mul_add(y, c);
        let mut step = 0.0;
        if f.abs() < df.abs() * 0.1 {
            step = f / df;
        }
        x[i] -= step;
    }
    n
}

/// Find the real zeros of a cubic polynomial p(x) = ax^3 + bx^2 + cx + d
pub(super) fn fpcuro(a: f64, b: f64, c: f64, d: f64, x: ArrayViewMut1<f64>) -> usize {
    let a1 = a.abs();
    let b1 = b.abs();
    let c1 = c.abs();
    let d1 = d.abs();
    if b1.max(c1).max(d1) < a1 * 1e4 {
        // 3rd degree polynomial
        real_zeros_cubic(a, b, c, d, x)
    } else if c1.max(d1) < b1 * 1e4 {
        // 2nd degree polynomial
        real_zeros_quadratic(a, b, c, d, x)
    } else if d1 < c1 * 1e4 {
        // 1st degree polynomial
        real_zeros_linear(a, b, c, d, x)
    } else {
        // constant function
        0
    }
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
    #[test]
    fn test_cubic() {
        let (a, b, c, d) = (1.0f64, -6.0, 11.0, -6.0);
        let mut x = array![0.0, 0.0, 0.0, 0.0];
        let n = fpcuro(a, b, c, d, x.view_mut());
        println!("{:?}", x);
        assert!((x[0] - 1.0).abs() < 1e-5);
        assert!((x[2] - 2.0).abs() < 1e-5);
        assert!((x[1] - 3.0).abs() < 1e-5);
    }
}
