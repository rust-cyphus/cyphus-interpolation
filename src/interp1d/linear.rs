use super::acc::InterpAccel;
use super::util::{accel_find, bsearch};
use crate::traits::{Interp1d, Interp1dAcc};

/// Linear interpolator
pub struct LinearInterp {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
}

impl LinearInterp {
    /// Create a new linear interpolator for a given abscissa `x` and
    /// ordinates `y`.
    pub fn new(x: &[f64], y: &[f64]) -> Self {
        LinearInterp {
            x: x.to_owned(),
            y: y.to_owned(),
        }
    }
    fn bsearch(&self, x: f64) -> usize {
        bsearch(&self.x, x, 0, self.x.len() - 1)
    }
}

impl Interp1d for LinearInterp {
    /// Evaluate the spline at `x`.
    fn eval(&self, x: f64) -> f64 {
        if x < self.x[0] || x > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        let idx = self.bsearch(x);

        let x_l = self.x[idx];
        let x_h = self.x[idx + 1];
        let y_l = self.y[idx];
        let y_h = self.y[idx + 1];
        let dx = x_h - x_l;

        if dx > 0.0 {
            y_l + (x - x_l) / dx * (y_h - y_l)
        } else {
            f64::NAN
        }
    }
    /// Evaluate the first derivative of the spline at `x`.
    fn derivative(&self, n: usize, x: f64) -> f64 {
        match n {
            0 => self.eval(x),
            1 => {
                if x < self.x[0] || x > self.x[self.x.len() - 1] {
                    return f64::NAN;
                }
                let idx = self.bsearch(x);

                let x_l = self.x[idx];
                let x_h = self.x[idx + 1];
                let y_l = self.y[idx];
                let y_h = self.y[idx + 1];
                let dx = x_h - x_l;
                let dy = y_h - y_l;

                if dx > 0.0 {
                    dy / dx
                } else {
                    f64::NAN
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
        if a > b || a < self.x[0] || b > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        let idx_a = self.bsearch(a);
        let idx_b = self.bsearch(b);

        let mut res = 0.0;
        for i in idx_a..(idx_b + 1) {
            let x_l = self.x[i];
            let x_h = self.x[i + 1];
            let y_l = self.y[i];
            let y_h = self.y[i + 1];
            let dx = x_h - x_l;

            if dx != 0.0 {
                if i == idx_a || i == idx_b {
                    let x1 = if i == idx_a { a } else { x_l };
                    let x2 = if i == idx_b { b } else { x_h };
                    let d = (y_h - y_l) / dx;
                    res += (x2 - x1) * (y_l + 0.5 * d * ((x2 - x_l) + (x1 - x_l)));
                } else {
                    res += 0.5 * dx * (y_l + y_h);
                }
            }
        }
        res
    }
}

/// Linear interpolator which uses an internal accelerator for quicker
/// evaluation. Requires exclusive reference.
pub struct LinearInterpAcc {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
    /// Accelerator
    acc: InterpAccel,
}

impl LinearInterpAcc {
    /// Create a new linear interpolator for a given abscissa `x` and
    /// ordinates `y`.
    pub fn new(x: &[f64], y: &[f64]) -> Self {
        Self {
            x: x.to_owned(),
            y: y.to_owned(),
            acc: InterpAccel::new(),
        }
    }
    fn bsearch(&self, x: f64) -> usize {
        bsearch(&self.x, x, 0, self.x.len() - 1)
    }
    fn accel_find(&mut self, x: f64) -> usize {
        accel_find(&self.x, x, &mut self.acc)
    }
}

impl Interp1dAcc for LinearInterpAcc {
    /// Evaluate the spline at `x` using the internal accelerator.
    fn eval(&mut self, x: f64) -> f64 {
        if x < self.x[0] || x > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        let idx = self.accel_find(x);

        let x_l = self.x[idx];
        let x_h = self.x[idx + 1];
        let y_l = self.y[idx];
        let y_h = self.y[idx + 1];
        let dx = x_h - x_l;

        if dx > 0.0 {
            y_l + (x - x_l) / dx * (y_h - y_l)
        } else {
            f64::NAN
        }
    }
    /// Evaluate the first derivative of the spline at `x` using the internal
    /// accelerator.
    fn derivative(&mut self, n: usize, x: f64) -> f64 {
        match n {
            0 => self.eval(x),
            1 => {
                if x < self.x[0] || x > self.x[self.x.len() - 1] {
                    return f64::NAN;
                }

                let idx = self.accel_find(x);

                let x_l = self.x[idx];
                let x_h = self.x[idx + 1];
                let y_l = self.y[idx];
                let y_h = self.y[idx + 1];
                let dx = x_h - x_l;
                let dy = y_h - y_l;

                if dx > 0.0 {
                    dy / dx
                } else {
                    f64::NAN
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
        if a > b || a < self.x[0] || b > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        let idx_a = self.accel_find(a);
        let idx_b = self.accel_find(b);

        let mut res = 0.0;
        for i in idx_a..(idx_b + 1) {
            let x_l = self.x[i];
            let x_h = self.x[i + 1];
            let y_l = self.y[i];
            let y_h = self.y[i + 1];
            let dx = x_h - x_l;

            if dx != 0.0 {
                if i == idx_a || i == idx_b {
                    let x1 = if i == idx_a { a } else { x_l };
                    let x2 = if i == idx_b { b } else { x_h };
                    let d = (y_h - y_l) / dx;
                    res += (x2 - x1) * (y_l + 0.5 * d * ((x2 - x_l) + (x1 - x_l)));
                } else {
                    res += 0.5 * dx * (y_l + y_h);
                }
            }
        }
        res
    }
}
