use super::acc::InterpAccel;
use super::util::{accel_find, bsearch, integ_eval, solve_tridiag};
use crate::traits::{Interp1d, Interp1dAcc};

/// Cubic spline interpolator
pub struct CubicSpline {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
    /// Coefficients
    c: Vec<f64>,
}

/// Cubic spline interpolator which uses an internal accelerator to speed up
/// evaluation. Requires exclusive reference for evaluations.
pub struct CubicSplineAcc {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
    /// Coefficients
    c: Vec<f64>,
    /// Accelerator
    acc: InterpAccel,
}

impl CubicSpline {
    /// Create a new spline interpolator using cubic interpolants from
    /// abscissa `x` and ordinates `y`.
    pub fn new(x: &[f64], y: &[f64]) -> Self {
        let npts = x.len();
        let max_idx = npts - 1;
        let sys_size = max_idx - 1;

        let mut c = vec![0.0; npts];
        let mut gg = vec![0.0; npts];
        let mut diag = vec![0.0; npts];
        let mut offdiag = vec![0.0; npts];

        for i in 0..sys_size {
            let h = x[i + 1] - x[i];
            let hp1 = x[i + 2] - x[i + 1];
            let ydiff = y[i + 1] - y[i];
            let ydiffpt = y[i + 2] - y[i + 1];
            let g = if h != 0.0 { 1.0 / h } else { 0.0 };
            let gp1 = if hp1 != 0.0 { 1.0 / hp1 } else { 0.0 };

            offdiag[i] = hp1;
            diag[i] = 2.0 * (hp1 + h);
            gg[i] = 3.0 * (ydiffpt * gp1 - ydiff * g);
        }

        if sys_size == 1 {
            c[1] = gg[0] / diag[0];
        } else {
            solve_tridiag(&diag, &offdiag, &gg, &mut c[1..], sys_size);
        }
        Self {
            x: x.to_owned(),
            y: y.to_owned(),
            c,
        }
    }
    fn coeff_calc(&self, dy: f64, dx: f64, idx: usize) -> (f64, f64, f64) {
        let c = self.c[idx];
        let cp1 = self.c[idx + 1];
        let b = dy / dx - dx * (cp1 + 2.0 * c) / 3.0;
        let d = (cp1 - c) / (3.0 * dx);

        (b, c, d)
    }
    fn bsearch(&self, x: f64) -> usize {
        bsearch(&self.x, x, 0, self.x.len() - 1)
    }
}

impl Interp1d for CubicSpline {
    /// Evaluate the spline at the point `x`.
    fn eval(&self, x: f64) -> f64 {
        let idx = self.bsearch(x);

        let xh = self.x[idx + 1];
        let xl = self.x[idx];
        let dx = xh - xl;

        if dx > 0.0 {
            let yh = self.y[idx + 1];
            let yl = self.y[idx];
            let dy = yh - yl;
            let delx = x - xl;
            let (b, c, d) = self.coeff_calc(dy, dx, idx);
            d.mul_add(delx, c).mul_add(delx, b).mul_add(delx, yl)
        } else {
            0.0
        }
    }
    /// Evaluate the first derivative of the spline at the point `x`.
    fn derivative(&self, n: usize, x: f64) -> f64 {
        match n {
            0 => self.eval(x),
            1 => {
                let idx = self.bsearch(x);

                let xh = self.x[idx + 1];
                let xl = self.x[idx];
                let dx = xh - xl;

                if dx > 0.0 {
                    let yh = self.y[idx + 1];
                    let yl = self.y[idx];
                    let dy = yh - yl;
                    let delx = x - xl;
                    let (b, c, d) = self.coeff_calc(dy, dx, idx);
                    (3.0 * d).mul_add(delx, 2.0 * c).mul_add(delx, b)
                } else {
                    0.0
                }
            }
            2 => {
                let idx = self.bsearch(x);

                let xh = self.x[idx + 1];
                let xl = self.x[idx];
                let dx = xh - xl;

                if dx > 0.0 {
                    let yl = self.y[idx];
                    let yh = self.y[idx + 1];
                    let dy = yh - yl;
                    let delx = x - xl;
                    let (_, c, d) = self.coeff_calc(dy, dx, idx);
                    delx.mul_add(6.0 * d, 2.0 * c)
                } else {
                    0.0
                }
            }
            _ => {
                if x < self.x[0] || x > self.x[self.x.len() - 1] {
                    return f64::NAN;
                }
                0.0
            }
        }
    }
    /// Integrate the spline from `a` to `b`.
    fn integrate(&self, a: f64, b: f64) -> f64 {
        let idx_a = self.bsearch(a);
        let idx_b = self.bsearch(b);

        let mut result = 0.0;
        for i in idx_a..(idx_b + 1) {
            let xh = self.x[i + 1];
            let xl = self.x[i];
            let yh = self.y[i + 1];
            let yl = self.y[i];
            let dx = xh - xl;
            let dy = yh - yl;
            if dx != 0.0 {
                let (bi, ci, di) = self.coeff_calc(dy, dx, i);
                result += if i == idx_a || i == idx_b {
                    let x1 = if i == idx_a { a } else { xl };
                    let x2 = if i == idx_b { b } else { xh };
                    integ_eval(yl, bi, ci, di, xl, x1, x2)
                } else {
                    (0.25 * di)
                        .mul_add(dx, ci / 3.0)
                        .mul_add(dx, 0.5 * bi)
                        .mul_add(dx, yl)
                        * dx
                };
            } else {
                return 0.0;
            }
        }
        result
    }
}

impl CubicSplineAcc {
    /// Create a new spline interpolator using cubic interpolants from
    /// abscissa `x` and ordinates `y`.
    pub fn new(x: &[f64], y: &[f64]) -> Self {
        let npts = x.len();
        let max_idx = npts - 1;
        let sys_size = max_idx - 1;

        let mut c = vec![0.0; npts];
        let mut gg = vec![0.0; npts];
        let mut diag = vec![0.0; npts];
        let mut offdiag = vec![0.0; npts];

        for i in 0..sys_size {
            let h = x[i + 1] - x[i];
            let hp1 = x[i + 2] - x[i + 1];
            let ydiff = y[i + 1] - y[i];
            let ydiffpt = y[i + 2] - y[i + 1];
            let g = if h != 0.0 { 1.0 / h } else { 0.0 };
            let gp1 = if hp1 != 0.0 { 1.0 / hp1 } else { 0.0 };

            offdiag[i] = hp1;
            diag[i] = 2.0 * (hp1 + h);
            gg[i] = 3.0 * (ydiffpt * gp1 - ydiff * g);
        }

        if sys_size == 1 {
            c[1] = gg[0] / diag[0];
        } else {
            solve_tridiag(&diag, &offdiag, &gg, &mut c[1..], sys_size);
        }
        Self {
            x: x.to_owned(),
            y: y.to_owned(),
            c,
            acc: InterpAccel::new(),
        }
    }
    fn coeff_calc(&self, dy: f64, dx: f64, idx: usize) -> (f64, f64, f64) {
        let c = self.c[idx];
        let cp1 = self.c[idx + 1];
        let b = dy / dx - dx * (cp1 + 2.0 * c) / 3.0;
        let d = (cp1 - c) / (3.0 * dx);

        (b, c, d)
    }
    fn bsearch(&self, x: f64) -> usize {
        bsearch(&self.x, x, 0, self.x.len() - 1)
    }
    fn accel_find(&mut self, x: f64) -> usize {
        accel_find(&self.x, x, &mut self.acc)
    }
}

impl Interp1dAcc for CubicSplineAcc {
    /// Evaluate the spline at the point `x` using the internal accelerator.
    fn eval(&mut self, x: f64) -> f64 {
        let idx = self.accel_find(x);

        let xh = self.x[idx + 1];
        let xl = self.x[idx];
        let dx = xh - xl;

        if dx > 0.0 {
            let yh = self.y[idx + 1];
            let yl = self.y[idx];
            let dy = yh - yl;
            let delx = x - xl;
            let (b, c, d) = self.coeff_calc(dy, dx, idx);
            d.mul_add(delx, c).mul_add(delx, b).mul_add(delx, yl)
        } else {
            0.0
        }
    }
    /// Evaluate the derivative of the spline at the point `x` using an
    /// internal accelerator.
    fn derivative(&mut self, n: usize, x: f64) -> f64 {
        match n {
            0 => self.eval(x),
            1 => {
                let idx = self.accel_find(x);

                let xh = self.x[idx + 1];
                let xl = self.x[idx];
                let dx = xh - xl;

                if dx > 0.0 {
                    let yh = self.y[idx + 1];
                    let yl = self.y[idx];
                    let dy = yh - yl;
                    let delx = x - xl;
                    let (b, c, d) = self.coeff_calc(dy, dx, idx);
                    (3.0 * d).mul_add(delx, 2.0 * c).mul_add(delx, b)
                } else {
                    0.0
                }
            }
            2 => {
                let idx = self.accel_find(x);

                let xh = self.x[idx + 1];
                let xl = self.x[idx];
                let dx = xh - xl;

                if dx > 0.0 {
                    let yl = self.y[idx];
                    let yh = self.y[idx + 1];
                    let dy = yh - yl;
                    let delx = x - xl;
                    let (_, c, d) = self.coeff_calc(dy, dx, idx);
                    delx.mul_add(6.0 * d, 2.0 * c)
                } else {
                    0.0
                }
            }
            _ => {
                if x < self.x[0] || x > self.x[self.x.len() - 1] {
                    return f64::NAN;
                }
                0.0
            }
        }
    }
    /// Integrate the spline from `a` to `b` using the internal accelerator.
    fn integrate(&mut self, a: f64, b: f64) -> f64 {
        let idx_a = self.accel_find(a);
        let idx_b = self.accel_find(b);

        let mut result = 0.0;
        for i in idx_a..(idx_b + 1) {
            let xh = self.x[i + 1];
            let xl = self.x[i];
            let yh = self.y[i + 1];
            let yl = self.y[i];
            let dx = xh - xl;
            let dy = yh - yl;
            if dx != 0.0 {
                let (bi, ci, di) = self.coeff_calc(dy, dx, i);
                result += if i == idx_a || i == idx_b {
                    let x1 = if i == idx_a { a } else { xl };
                    let x2 = if i == idx_b { b } else { xh };
                    integ_eval(yl, bi, ci, di, xl, x1, x2)
                } else {
                    (0.25 * di)
                        .mul_add(dx, ci / 3.0)
                        .mul_add(dx, 0.5 * bi)
                        .mul_add(dx, yl)
                        * dx
                };
            } else {
                return 0.0;
            }
        }
        result
    }
}
