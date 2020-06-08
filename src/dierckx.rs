//pub(super) mod curfit;
pub(super) mod fpader;
pub(super) mod fpback;
pub(super) mod fpbspl;
pub(super) mod fpchec;
//pub(super) mod fpcurf;
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

pub struct UnivariateSpline {
    /// 1-D array of indpendent input data. Must be increasing; must be
    /// strictly increasing if `s` is 0.
    pub x: Array1<f64>,
    /// 1-D array of dependent input data, of the same length as `x`.
    pub y: Array1<f64>,
    /// Weights for spline fitting. Must be positvie. If none (default),
    /// weights are all equal.
    pub w: Array1<f64>,
    /// Tuple specifying the boundary of the approximation interval. If None
    /// (default), `bbox=[x[0], x[x.len()-1]]`
    pub bbox: (f64, f64),
    /// Degree of the smoothing spline. Must be <= 5. Default is `k=3`, a
    /// cubic spline.
    pub k: usize,
    /// Positive smoothing factor used to choose the number of knots. Number
    /// of knots will be increased until the smoothing condition is
    /// satisfied: sum((w[i] * (y[i] - spl(x[i]))) <= s. In None (default),
    /// `s = w.len()` which should be a good valud if `1/w[i]` is an
    /// estimate of the standard deviation of `y[i]`. If 0, spline will
    /// interpolate through all data points.
    pub s: f64,
    /// Controls the extrapolation mode for elements not in the interval
    /// defined by the knot sequence.
    pub ext: usize,
    // Number of knots.
    n: usize,
    // Size of x, y and w
    m: usize,
    // Size of t, c fpint and z
    nest: usize,
    // Width of a and q
    k1: usize,
    // Width of b and g
    k2: usize,
    // Knot locations
    t: Array1<f64>,
    // Spline coefficients
    c: Array1<f64>,
    // Working arrays for curfit
    fpint: Array1<f64>,
    z: Array1<f64>,
    a: Array2<f64>,
    b: Array2<f64>,
    g: Array2<f64>,
    q: Array2<f64>,
    iwrk: Array1<i32>,
}
