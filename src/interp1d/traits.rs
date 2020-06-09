pub trait Interp1d {
    /// Evaluate the interpolator at a point
    fn eval(&self, x: f64) -> f64;
    /// Evaluate interpolator at a point using accelerator. Requires
    /// interpolator to be mutable.
    fn eval_acc(&mut self, x: f64) -> f64;
    /// Evalaute the derivative of the interpolator at a point.
    fn deriv(&self, x: f64) -> f64;
    /// Evalaute the derivative of the interpolator at a point using
    /// accelerator. Requires interpolator to be mutable.
    fn deriv_acc(&mut self, x: f64) -> f64;
    /// Evalaute the 2nd derivative of the interpolator at a point.
    fn second_deriv(&self, x: f64) -> f64;
    /// Evalaute the 2nd derivative of the interpolator at a point using
    /// accelerator. Requires interpolator to be mutable.
    fn second_deriv_acc(&mut self, x: f64) -> f64;
    /// Compute the integral of an interpolator from `a` to `b`.
    fn integrate(&self, a: f64, b: f64) -> f64;
    /// Compute the integral of an interpolator from `a` to `b` using
    /// accelerator. Requires interpolator to be mutable.    
    fn integrate_acc(&mut self, a: f64, b: f64) -> f64;
}
