pub mod acc;
pub mod cubic;
pub mod linear;
pub(crate) mod traits;
pub(crate) mod util;

pub use acc::InterpAccel;
pub use cubic::*;
pub use linear::*;
pub use traits::Interp1d;

/// Cubic spline interpolator
pub struct CubicSpline {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
    /// Coefficients
    c: Vec<f64>,
    /// Accelerator
    acc: InterpAccel,
}

/// Linear interpolator
pub struct LinearInterp {
    /// x-data
    x: Vec<f64>,
    /// y-data
    y: Vec<f64>,
    /// Accelerator
    acc: InterpAccel,
}
