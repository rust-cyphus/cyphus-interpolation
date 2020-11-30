use cyphus_interpolation::prelude::*;
use ndarray::prelude::*;

#[test]
fn test_eval() {
    let npts = 50;
    let interval = (0.0, std::f64::consts::PI);
    let xs = Array1::linspace(interval.0, interval.1, npts);
    let ys = xs.mapv(f64::sin);

    let spline = UnivariateSplineBuilder::default(&xs, &ys).build().unwrap();
    let ys = spline.eval_array(&xs);
    for (x, y) in xs.iter().zip(ys.iter()) {
        assert!((*y - (*x).sin()).abs() < 1e-10);
    }
}
#[test]
fn test_extrapolation_0() {
    let npts = 50;
    let interval = (0.0, std::f64::consts::PI);
    let xs = Array1::linspace(interval.0, interval.1, npts);
    let ys = xs.mapv(f64::cos);

    let spline = UnivariateSplineBuilder::default(&xs, &ys)
        .extrapolation(0)
        .build()
        .unwrap();
    dbg!(spline.eval(-1.0));
}
#[test]
/// Test that when ext = 1 out-of-bounds evaluation returns zero
fn test_extrapolation_1() {
    let npts = 50;
    let interval = (0.0, std::f64::consts::PI);
    let xs = Array1::linspace(interval.0, interval.1, npts);
    let ys = xs.mapv(f64::cos);

    let spline = UnivariateSplineBuilder::default(&xs, &ys)
        .extrapolation(1)
        .build()
        .unwrap();
    assert!(spline.eval(xs[0] - 1.0) == 0.0);
    assert!(spline.eval(xs[npts - 1] + 1.0) == 0.0);
}
#[test]
#[should_panic]
/// Test that when ext = 2 out-of-bounds evaluation throws panics    
fn test_extrapolation_2() {
    let npts = 50;
    let interval = (0.0, std::f64::consts::PI);
    let xs = Array1::linspace(interval.0, interval.1, npts);
    let ys = xs.mapv(f64::cos);

    let spline = UnivariateSplineBuilder::default(&xs, &ys)
        .extrapolation(2)
        .build()
        .unwrap();
    spline.eval(-1.0);
}
#[test]
/// Test that when ext = 3 out-of-bounds evaluation returns boundary
/// value
fn test_extrapolation_3() {
    let npts = 50;
    let interval = (0.0, std::f64::consts::PI);
    let xs = Array1::linspace(interval.0, interval.1, npts);
    let ys = xs.mapv(f64::cos);

    let spline = UnivariateSplineBuilder::default(&xs, &ys)
        .extrapolation(3)
        .build()
        .unwrap();
    assert!((spline.eval(interval.0 - 1.0) - ys[0]).abs() < 1e-10);
    assert!((spline.eval(interval.1 + 1.0) - ys[npts - 1]).abs() < 1e-10);
}
#[test]
fn test_integrate() {
    let npts = 100;
    let interval = (0.0, std::f64::consts::PI);
    let xs = Array1::linspace(interval.0, interval.1, npts);
    let ys = xs.mapv(f64::sin);

    let spline = UnivariateSplineBuilder::default(&xs, &ys).build().unwrap();
    let int = spline.integrate(1.0, 2.0);
    assert!((int - 0.956449142415282).abs() < 1e-5);
}
#[test]
fn test_derivatives() {
    let npts = 100;
    let interval = (0.0, std::f64::consts::PI);
    let step = (interval.1 - interval.0) / (npts - 1) as f64;
    let mut xs = Array1::<f64>::zeros(npts);
    let mut ys = Array1::<f64>::zeros(npts);

    for i in 0..npts {
        xs[i] = interval.0 + step * i as f64;
        ys[i] = xs[i].sin();
    }

    let spline = UnivariateSplineBuilder::default(&xs, &ys).build().unwrap();

    let dys = spline.derivative_array(1, &xs);
    for (x, dy) in xs.iter().zip(dys.iter()) {
        // let dy0 = (*x).cos();
        // println!("{}, {}, {}, {:e}", *x, *dy, dy0, (*dy - dy0).abs());
        assert!((*dy - (*x).cos()).abs() < 2e-7);
    }
}
#[test]
fn test_roots() {
    let pi = std::f64::consts::PI;
    let xs = Array::linspace(-3.0 * pi / 2.0, 3.0 * pi / 2.0, 100);
    let ys = xs.mapv(f64::sin);

    let spline = UnivariateSplineBuilder::default(&xs, &ys).build().unwrap();
    let zeros = spline.roots();

    assert!((zeros[0] + pi).abs() < 1e-5);
    assert!((zeros[1]).abs() < 1e-5);
    assert!((zeros[2] - pi).abs() < 1e-5);
}
