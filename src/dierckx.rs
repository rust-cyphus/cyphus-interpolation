//! # `dierckx.rs`
//! Port of the DIERCKX fortran library in the same style as `scipy`.
//!
//! ## Examples
//! [`UnivariateSpline`]
//!

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
use crate::traits::Interp1d;
use ndarray::prelude::*;

#[derive(Debug)]
pub enum UnivariateSplineBuildErr {
    StorageErr,
    ImpossibleResultErr,
    MaxIterErr,
    InvalidDataErr,
}

#[derive(Debug, Clone)]
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
    /// Create a `UnivariateSplineBuilder` with abscissa `x` and ordinates `y`.
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
    /// Set the weights of the abscissa.
    pub fn weights(mut self, w: Array1<f64>) -> UnivariateSplineBuilder {
        self.w = Some(w);
        self
    }
    /// Set the boundary of the approximation interval.
    pub fn boundary(mut self, bbox: (f64, f64)) -> UnivariateSplineBuilder {
        self.bbox = Some(bbox);
        self
    }
    /// Set the degree of the spline.
    pub fn degree(mut self, k: usize) -> UnivariateSplineBuilder {
        self.k = Some(k);
        self
    }
    /// Set the smoothing factor of the spline.
    pub fn smoothing_factor(mut self, s: f64) -> UnivariateSplineBuilder {
        self.s = Some(s);
        self
    }
    /// Specify the extrapolation type:
    /// - if ext=0, return the extrapolated value,
    /// - if ext=1, return 0,
    /// - if ext=2, raise an error,
    /// - if ext=3, return the boundary value.
    pub fn extrapolation(mut self, ext: usize) -> UnivariateSplineBuilder {
        self.ext = Some(ext);
        self
    }
    pub fn build(self) -> Result<UnivariateSpline, UnivariateSplineBuildErr> {
        let m = self.x.len();
        let w = self.w.unwrap_or_else(|| Array1::<f64>::ones(m));
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

impl Interp1d for UnivariateSpline {
    /// Evaluate the spline at a single point.
    fn eval(&self, x: f64) -> f64 {
        let mut ier = 0;
        let y = splev::splev_point(
            self.t.view(),
            self.n,
            self.c.view(),
            self.k,
            x,
            self.ext,
            &mut ier,
        );
        match ier {
            0 => y,
            1 => panic!("Extrapolation out of bounds."),
            10 => panic!("Invalid data."),
            _ => unreachable!(),
        }
    }
    fn derivative(&self, nu: usize, x: f64) -> f64 {
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
    fn integrate(&self, a: f64, b: f64) -> f64 {
        splint::splint(self.t.view(), self.n, self.c.view(), self.k, a, b)
    }
}

impl UnivariateSpline {
    /// Evaluate the spline at a vector of points.
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
        match ier {
            0 => y,
            1 => panic!("Extrapolation out of bounds."),
            10 => panic!("Invalid data."),
            _ => unreachable!(),
        }
    }
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
    pub fn roots(&self) -> Array1<f64> {
        let mut ier = 0;
        let zero = sproot::sproot(self.t.view(), self.n, self.c.view(), &mut ier);
        dbg!(ier);
        zero
    }
}
