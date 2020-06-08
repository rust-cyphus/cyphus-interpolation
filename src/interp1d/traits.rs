pub trait Interp {
    fn eval(&mut self, x: f64) -> f64;
    fn deriv(&mut self, x: f64) -> f64;
    fn second_deriv(&mut self, x: f64) -> f64;
    fn integrate(&mut self, a: f64, b: f64) -> f64;
}
