pub(super) mod curfit;
pub(super) mod fpader;
pub(super) mod fpback;
pub(super) mod fpbspl;
pub(super) mod fpchec;
pub(super) mod fpcurf;
pub(super) mod fpcuro;
pub(super) mod fpdisc;
pub(super) mod fpgivs;
pub(super) mod fpintb;
pub(super) mod fpknot;
pub(super) mod fprati;
pub(super) mod fprota;
pub(super) mod spalde;
pub(super) mod splder;
pub(super) mod splev;
pub(super) mod splint;
pub(super) mod sproot;
use ndarray::prelude::*;

#[derive(Debug)]
pub enum UnivariateSplineBuildErr {
    StorageErr,
    ImpossibleResultErr,
    MaxIterErr,
    InvalidDataErr,
}

pub struct UnivariateSpline {
    // Degree of the smoothing spline. Must be <= 5. Default is `k=3`, a
    // cubic spline.
    pub(crate) k: usize,
    // Controls the extrapolation mode for elements not in the interval
    // defined by the knot sequence.
    pub(crate) ext: usize,
    // Number of knots.
    pub(crate) n: usize,
    // Size of x, y and w
    pub(crate) m: usize,
    // Size of t, c fpint and z
    pub(crate) nest: usize,
    // Width of a and q
    pub(crate) k1: usize,
    // Width of b and g
    pub(crate) k2: usize,
    // Knot locations
    pub(crate) t: Array1<f64>,
    // Spline coefficients
    pub(crate) c: Array1<f64>,
}

pub struct UnivariateSplineBuilder {
    /// 1-D array of indpendent input data. Must be increasing; must be
    /// strictly increasing if `s` is 0.
    pub x: Array1<f64>,
    /// 1-D array of dependent input data, of the same length as `x`.
    pub y: Array1<f64>,
    /// Weights for spline fitting. Must be positvie. If none (default),
    /// weights are all equal.
    pub w: Option<Array1<f64>>,
    /// Tuple specifying the boundary of the approximation interval. If None
    /// (default), `bbox=[x[0], x[x.len()-1]]`
    pub bbox: Option<(f64, f64)>,
    /// Degree of the smoothing spline. Must be <= 5. Default is `k=3`, a
    /// cubic spline.
    pub k: Option<usize>,
    /// Positive smoothing factor used to choose the number of knots. Number
    /// of knots will be increased until the smoothing condition is
    /// satisfied: sum((w[i] * (y[i] - spl(x[i]))) <= s. In None (default),
    /// `s = w.len()` which should be a good valud if `1/w[i]` is an
    /// estimate of the standard deviation of `y[i]`. If 0, spline will
    /// interpolate through all data points.
    pub s: Option<f64>,
    /// Controls the extrapolation mode for elements not in the interval
    /// defined by the knot sequence.
    pub ext: Option<usize>,
}

impl UnivariateSplineBuilder {
    #[allow(dead_code)]
    pub fn default(x: &Array1<f64>, y: &Array1<f64>) -> UnivariateSplineBuilder {
        UnivariateSplineBuilder {
            x: x.clone(),
            y: y.clone(),
            w: None,
            bbox: None,
            k: None,
            s: None,
            ext: None,
        }
    }
    #[allow(dead_code)]
    pub fn weights(&mut self, w: Array1<f64>) -> &mut UnivariateSplineBuilder {
        self.w = Some(w.clone());
        self
    }
    #[allow(dead_code)]
    pub fn boundary(&mut self, bbox: (f64, f64)) -> &mut UnivariateSplineBuilder {
        self.bbox = Some(bbox.clone());
        self
    }
    #[allow(dead_code)]
    pub fn degree(&mut self, k: usize) -> &mut UnivariateSplineBuilder {
        self.k = Some(k);
        self
    }
    #[allow(dead_code)]
    pub fn smoothing_factor(&mut self, s: f64) -> &mut UnivariateSplineBuilder {
        self.s = Some(s);
        self
    }
    #[allow(dead_code)]
    pub fn extrapolation(&mut self, ext: usize) -> &mut UnivariateSplineBuilder {
        self.ext = Some(ext);
        self
    }
    #[allow(dead_code)]
    pub fn build(self) -> Result<UnivariateSpline, UnivariateSplineBuildErr> {
        let m = self.x.len();
        let w = self.w.unwrap_or(Array1::<f64>::ones(m));
        let bbox = self.bbox.unwrap_or((self.x[0], self.x[m - 1]));
        let mut k = self.k.unwrap_or(3);
        let s = self.s.unwrap_or(0.0);
        let ext = self.ext.unwrap_or(0);

        let nest = m + k + 1;
        let k1 = k + 1;
        let k2 = k1 + 1;

        let mut n = 0;
        let mut fp = 0.0;
        let mut ier = 0;

        // Knot locations
        let mut t = Array1::<f64>::zeros(nest);
        // Spline coefficients
        let mut c = Array1::<f64>::zeros(nest);
        // Working arrays for curfit
        let mut fpint = Array1::<f64>::zeros(nest);
        let mut z = Array1::<f64>::zeros(nest);
        let mut a = Array2::<f64>::zeros((nest, k1));
        let mut b = Array2::<f64>::zeros((nest, k2));
        let mut g = Array2::<f64>::zeros((nest, k2));
        let mut q = Array2::<f64>::zeros((m, k1));
        let mut iwrk = Array1::<usize>::zeros(nest);

        let mut iopt = 0;

        curfit::curfit(
            &mut iopt,
            m,
            self.x.view(),
            self.y.view(),
            w.view(),
            bbox.0,
            bbox.1,
            &mut k,
            s,
            nest,
            &mut n,
            t.view_mut(),
            c.view_mut(),
            &mut fp,
            fpint.view_mut(),
            z.view_mut(),
            a.view_mut(),
            b.view_mut(),
            g.view_mut(),
            q.view_mut(),
            iwrk.view_mut(),
            &mut ier,
        );

        if ier > 0 {
            if ier == 1 {
                Err(UnivariateSplineBuildErr::StorageErr)
            } else if ier == 2 {
                Err(UnivariateSplineBuildErr::ImpossibleResultErr)
            } else if ier == 3 {
                Err(UnivariateSplineBuildErr::MaxIterErr)
            } else {
                Err(UnivariateSplineBuildErr::InvalidDataErr)
            }
        } else {
            Ok(UnivariateSpline {
                // Degree of interpolation
                k,
                // Extrapolation type
                ext,
                // Number of knots.
                n,
                // Size of x, y and w
                m,
                // Size of t, c fpint and z
                nest,
                // Width of a and q
                k1,
                // Width of b and g
                k2,
                // Knot locations
                t,
                // Spline coefficients
                c,
            })
        }
    }
}

impl UnivariateSpline {
    /// Evaluate the spline at a single point.
    #[allow(dead_code)]
    pub fn eval(&self, x: f64) -> f64 {
        let mut ier = 0;
        splev::splev_point(
            self.t.view(),
            self.n,
            self.c.view(),
            self.k,
            x,
            self.ext,
            &mut ier,
        )
    }
    /// Evaluate the spline at a vector of points.
    #[allow(dead_code)]
    pub fn eval_array(&self, x: &Array1<f64>) -> Array1<f64> {
        let mut ier = 0;
        let mut y = Array1::<f64>::zeros(x.len());
        let m = y.len();
        splev::splev(
            self.t.view(),
            self.n,
            self.c.view(),
            self.k,
            x.view(),
            y.view_mut(),
            m,
            self.ext,
            &mut ier,
        );
        y
    }
    #[allow(dead_code)]
    pub fn integrate(&self, a: f64, b: f64) -> f64 {
        splint::splint(self.t.view(), self.n, self.c.view(), self.k, a, b)
    }
    #[allow(dead_code)]
    pub fn derivative(&self, nu: usize, x: f64) -> f64 {
        if nu > self.k {}
        let mut ier = 0;
        let mut wrk = Array1::<f64>::zeros(self.n);
        let xx = array![x];
        let mut y = Array1::<f64>::zeros(1);

        splder::splder(
            self.t.view(),
            self.n,
            self.c.view(),
            self.k,
            nu,
            xx.view(),
            y.view_mut(),
            1,
            self.ext,
            wrk.view_mut(),
            &mut ier,
        );
        y[0]
    }
    #[allow(dead_code)]
    pub fn derivative_array(&self, nu: usize, x: &Array1<f64>) -> Array1<f64> {
        if nu > self.k {}
        let mut ier = 0;
        let mut wrk = Array1::<f64>::zeros(self.n);
        let mut y = Array1::<f64>::zeros(x.len());

        splder::splder(
            self.t.view(),
            self.n,
            self.c.view(),
            self.k,
            nu,
            x.view(),
            y.view_mut(),
            x.len(),
            self.ext,
            wrk.view_mut(),
            &mut ier,
        );
        y
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_eval() {
        let npts = 50;
        let interval = (0.0, std::f64::consts::PI);
        let step = (interval.1 - interval.0) / (npts - 1) as f64;
        let mut xs = Array1::<f64>::zeros(npts);
        let mut ys = Array1::<f64>::zeros(npts);

        for i in 0..npts {
            xs[i] = interval.0 + step * i as f64;
            ys[i] = xs[i].sin();
        }

        let spline = UnivariateSplineBuilder::default(&xs, &ys).build().unwrap();
        let ys = spline.eval_array(&xs);
        for (x, y) in xs.iter().zip(ys.iter()) {
            assert!((*y - (*x).sin()).abs() < 1e-10);
        }
    }
    #[test]
    fn test_integrate() {
        let npts = 50;
        let interval = (0.0, std::f64::consts::PI);
        let step = (interval.1 - interval.0) / (npts - 1) as f64;
        let mut xs = Array1::<f64>::zeros(npts);
        let mut ys = Array1::<f64>::zeros(npts);

        for i in 0..npts {
            xs[i] = interval.0 + step * i as f64;
            ys[i] = xs[i].sin();
        }

        let spline = UnivariateSplineBuilder::default(&xs, &ys).build().unwrap();
        let int = spline.integrate(1.0, 2.0);
        println!("{}", int);
        assert!((int - 0.956449142415282).abs() < 1e-5);
    }
    #[test]
    fn test_derivatives() {
        let npts = 50;
        let interval = (0.0, std::f64::consts::PI);
        let step = (interval.1 - interval.0) / (npts - 1) as f64;
        let mut xs = Array1::<f64>::zeros(npts);
        let mut ys = Array1::<f64>::zeros(npts);

        for i in 0..npts {
            xs[i] = interval.0 + step * i as f64;
            ys[i] = xs[i].sin();
        }

        let spline = UnivariateSplineBuilder::default(&xs, &ys).build().unwrap();

        println!("coeffs = {}", spline.c);
        println!("");
        println!("");
        println!("knots = {}", spline.t);
        println!("");
        println!("");

        let dys = spline.derivative_array(1, &xs);
        for (x, dy) in xs.iter().zip(dys.iter()) {
            let dy0 = (*x).cos();
            println!("{}, {}, {}, {}", *x, *dy, dy0, (*dy - dy0).abs());
            //assert!((*dy - (*x).cos()).abs() < 1e-10);
        }
        assert!(false);
    }
}
