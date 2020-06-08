pub(crate) mod acc;
pub mod cubic;
pub mod linear;
pub(crate) mod traits;
pub(crate) mod util;

/// Cubic spline interpolator
pub struct CubicSpline {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
    /// Coefficients
    c: Vec<f64>,
    /// Accelerator
    acc: acc::InterpAccel,
}

/// Linear interpolator
pub struct LinearInterp {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
    /// Accelerator
    acc: acc::InterpAccel,
}
