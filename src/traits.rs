/// Trait describing the functionality of a 1-dimensional interpolator.
pub trait Interp1d {
    /// Evaluate the interpolator at a point
    fn eval(&self, x: f64) -> f64;
    /// Evalaute the `n`th derivative of the interpolator at a point.
    fn derivative(&self, n: usize, x: f64) -> f64;
    /// Compute the integral of an interpolator from `a` to `b`.
    fn integrate(&self, a: f64, b: f64) -> f64;
}

/// Trait describing the functionality of a 1-dimensional interpolator which
/// uses an internal accelerator for improved speed.
pub trait Interp1dAcc {
    /// Evaluate interpolator at a point using accelerator. Requires
    /// interpolator to be mutable.
    fn eval(&mut self, x: f64) -> f64;
    /// Evalaute the `n`th derivative of the interpolator at a point using
    /// accelerator. Requires interpolator to be mutable.
    fn derivative(&mut self, n: usize, x: f64) -> f64;
    /// Compute the integral of an interpolator from `a` to `b` using
    /// accelerator. Requires interpolator to be mutable.    
    fn integrate(&mut self, a: f64, b: f64) -> f64;
}
