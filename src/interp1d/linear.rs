use super::acc::InterpAccel;
use super::traits::Interp;
use super::util::{accel_find, bsearch};
use super::LinearInterp;

impl LinearInterp {
    /// Create a new linear interpolator for a given set of `x` and `y` data.
    #[allow(dead_code)]
    pub fn new(x: &Vec<f64>, y: &Vec<f64>) -> LinearInterp {
        LinearInterp {
            x: x.clone(),
            y: y.clone(),
            acc: None,
        }
    }
    /// Create a new linear interpolator for a given set of `x` and `y` data.
    #[allow(dead_code)]
    pub fn with_acc(x: &Vec<f64>, y: &Vec<f64>) -> LinearInterp {
        LinearInterp {
            x: x.clone(),
            y: y.clone(),
            acc: Some(InterpAccel::new()),
        }
    }
    fn find(&mut self, x: f64) -> usize {
        match self.acc {
            Some(ref mut acc) => accel_find(&self.x, x, acc),
            None => bsearch(&self.x, x, 0, self.x.len() - 1),
        }
    }
}

impl Interp for LinearInterp {
    fn eval(&mut self, x: f64) -> f64 {
        if x < self.x[0] || x > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        let idx = self.find(x);

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
    fn deriv(&mut self, x: f64) -> f64 {
        if x < self.x[0] || x > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        let idx = self.find(x);

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
    fn second_deriv(&mut self, x: f64) -> f64 {
        if x < self.x[0] || x > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        0.0
    }
    fn integrate(&mut self, a: f64, b: f64) -> f64 {
        if a > b || a < self.x[0] || b > self.x[self.x.len() - 1] {
            return f64::NAN;
        }
        let idx_a = self.find(a);
        let idx_b = self.find(b);

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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_linear() {
        let data_x = vec![0.0, 1.0, 2.0, 3.0];
        let data_y = vec![0.0, 1.0, 2.0, 3.0];
        let test_x = vec![0.0, 0.5, 1.0, 1.5, 2.5, 3.0];
        let test_y = vec![0.0, 0.5, 1.0, 1.5, 2.5, 3.0];
        let test_dy = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let test_iy = vec![0.0, 0.125, 0.5, 9.0 / 8.0, 25.0 / 8.0, 9.0 / 2.0];

        let mut interp = LinearInterp::new(&data_x, &data_y);

        for i in 0..test_x.len() {
            let xt = test_x[i];
            let yt = test_y[i];
            let dyt = test_dy[i];
            let intt = test_iy[i];

            let y = interp.eval(xt);
            let dy = interp.deriv(xt);
            let int = interp.integrate(test_x[0], xt);

            let diff_y = y - yt;
            let diff_deriv = dy - dyt;
            let diff_int = int - intt;

            assert!(diff_y.abs() <= 1e-10);
            assert!(diff_deriv.abs() <= 1e-10);
            assert!(diff_int.abs() <= 1e-10);
        }
    }
}
