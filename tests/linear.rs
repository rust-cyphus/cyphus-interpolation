use cyphus_interpolation::prelude::*;

#[test]
fn test_linear() {
    let data_x = vec![0.0, 1.0, 2.0, 3.0];
    let data_y = vec![0.0, 1.0, 2.0, 3.0];
    let test_x = vec![0.0, 0.5, 1.0, 1.5, 2.5, 3.0];
    let test_y = vec![0.0, 0.5, 1.0, 1.5, 2.5, 3.0];
    let test_dy = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let test_iy = vec![0.0, 0.125, 0.5, 9.0 / 8.0, 25.0 / 8.0, 9.0 / 2.0];

    let interp = LinearInterp::new(&data_x, &data_y);

    for i in 0..test_x.len() {
        let xt = test_x[i];
        let yt = test_y[i];
        let dyt = test_dy[i];
        let intt = test_iy[i];

        let y = interp.eval(xt);
        let dy = interp.derivative(1, xt);
        let int = interp.integrate(test_x[0], xt);

        let diff_y = y - yt;
        let diff_deriv = dy - dyt;
        let diff_int = int - intt;

        assert!(diff_y.abs() <= 1e-10);
        assert!(diff_deriv.abs() <= 1e-10);
        assert!(diff_int.abs() <= 1e-10);
    }
}
