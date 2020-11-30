#![feature(test)]

extern crate test;
use cyphus_interpolation::prelude::*;
use ndarray::prelude::*;
use std::f64::consts::PI;
use test::black_box;

#[bench]
fn linear_eval(bench: &mut test::Bencher) {
    let npts = 20;
    let xs = Array1::linspace(0.0, PI, npts);
    let xs2 = Array1::linspace(0.0, PI, npts * 3);
    let ys = xs.mapv(f64::cos);
    let interp = LinearInterp::new(xs.as_slice().unwrap(), ys.as_slice().unwrap());

    let interp = black_box(interp);
    bench.iter(|| {
        for x in xs2.iter() {
            interp.eval(*x);
        }
    })
}

#[bench]
fn linear_acc_eval(bench: &mut test::Bencher) {
    let npts = 20;
    let xs = Array1::linspace(0.0, PI, npts);
    let xs2 = Array1::linspace(0.0, PI, npts * 3);
    let ys = xs.mapv(f64::cos);
    let interp = LinearInterpAcc::new(xs.as_slice().unwrap(), ys.as_slice().unwrap());

    let mut interp = black_box(interp);
    bench.iter(|| {
        for x in xs2.iter() {
            interp.eval(*x);
        }
    })
}

#[bench]
fn linear_univariate_spline_eval(bench: &mut test::Bencher) {
    let npts = 20;
    let xs = Array1::linspace(0.0, PI, npts);
    let xs2 = Array1::linspace(0.0, PI, npts * 3);
    let ys = xs.mapv(f64::cos);
    let interp = UnivariateSplineBuilder::default(&xs, &ys)
        .degree(1)
        .build()
        .unwrap();

    let interp = black_box(interp);
    bench.iter(|| {
        for x in xs2.iter() {
            interp.eval(*x);
        }
    })
}
