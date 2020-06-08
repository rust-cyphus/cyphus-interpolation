use super::acc::InterpAccel;
use super::traits::Interp;
use super::util::{accel_find, integ_eval, solve_tridiag};

pub struct CubicSpline {
    x: Vec<f64>,
    y: Vec<f64>,
    c: Vec<f64>,
    acc: InterpAccel,
}

impl CubicSpline {
    #[allow(dead_code)]
    pub fn new(x: &Vec<f64>, y: &Vec<f64>) -> CubicSpline {
        let npts = x.len();
        let max_idx = npts - 1;
        let sys_size = max_idx - 1;

        let mut c = vec![0.0; npts];
        let mut gg = vec![0.0; npts];
        let mut diag = vec![0.0; npts];
        let mut offdiag = vec![0.0; npts];

        for i in 0..sys_size {
            let h = x[i + 1] - x[i];
            let hp1 = x[i + 2] - x[i + 1];
            let ydiff = y[i + 1] - y[i];
            let ydiffpt = y[i + 2] - y[i + 1];
            let g = if h != 0.0 { 1.0 / h } else { 0.0 };
            let gp1 = if hp1 != 0.0 { 1.0 / hp1 } else { 0.0 };

            offdiag[i] = hp1;
            diag[i] = 2.0 * (hp1 + h);
            gg[i] = 3.0 * (ydiffpt * gp1 - ydiff * g);
        }

        if sys_size == 1 {
            c[1] = gg[0] / diag[0];
        } else {
            solve_tridiag(&diag, &offdiag, &gg, &mut c[1..], sys_size);
        }
        CubicSpline {
            x: x.clone(),
            y: y.clone(),
            c,
            acc: InterpAccel::new(),
        }
    }
    #[allow(dead_code)]
    fn coeff_calc(&self, dy: f64, dx: f64, idx: usize) -> (f64, f64, f64) {
        let c = self.c[idx];
        let cp1 = self.c[idx + 1];
        let b = dy / dx - dx * (cp1 + 2.0 * c) / 3.0;
        let d = (cp1 - c) / (3.0 * dx);

        (b, c, d)
    }
    #[allow(dead_code)]
    fn accel_find(&mut self, x: f64) -> usize {
        accel_find(&self.x, x, &mut self.acc)
    }
}

impl Interp for CubicSpline {
    fn eval(&mut self, x: f64) -> f64 {
        let idx = self.accel_find(x);

        let xh = self.x[idx + 1];
        let xl = self.x[idx];
        let dx = xh - xl;

        if dx > 0.0 {
            let yh = self.y[idx + 1];
            let yl = self.y[idx];
            let dy = yh - yl;
            let delx = x - xl;
            let (b, c, d) = self.coeff_calc(dy, dx, idx);
            d.mul_add(delx, c).mul_add(delx, b).mul_add(delx, yl)
        } else {
            0.0
        }
    }
    fn deriv(&mut self, x: f64) -> f64 {
        let idx = self.accel_find(x);

        let xh = self.x[idx + 1];
        let xl = self.x[idx];
        let dx = xh - xl;

        if dx > 0.0 {
            let yh = self.y[idx + 1];
            let yl = self.y[idx];
            let dy = yh - yl;
            let delx = x - xl;
            let (b, c, d) = self.coeff_calc(dy, dx, idx);
            (3.0 * d).mul_add(delx, 2.0 * c).mul_add(delx, b)
        } else {
            0.0
        }
    }
    fn second_deriv(&mut self, x: f64) -> f64 {
        let idx = self.accel_find(x);

        let xh = self.x[idx + 1];
        let xl = self.x[idx];
        let dx = xh - xl;

        if dx > 0.0 {
            let yl = self.y[idx];
            let yh = self.y[idx + 1];
            let dy = yh - yl;
            let delx = x - xl;
            let (_, c, d) = self.coeff_calc(dy, dx, idx);
            delx.mul_add(6.0 * d, 2.0 * c)
        } else {
            0.0
        }
    }
    fn integrate(&mut self, a: f64, b: f64) -> f64 {
        let idx_a = self.accel_find(a);
        let idx_b = self.accel_find(b);

        let mut result = 0.0;
        for i in idx_a..(idx_b + 1) {
            let xh = self.x[i + 1];
            let xl = self.x[i];
            let yh = self.y[i + 1];
            let yl = self.y[i];
            let dx = xh - xl;
            let dy = yh - yl;
            if dx != 0.0 {
                let (bi, ci, di) = self.coeff_calc(dy, dx, i);
                result += if i == idx_a || i == idx_b {
                    let x1 = if i == idx_a { a } else { xl };
                    let x2 = if i == idx_b { b } else { xh };
                    integ_eval(yl, bi, ci, di, xl, x1, x2)
                } else {
                    (0.25 * di)
                        .mul_add(dx, ci / 3.0)
                        .mul_add(dx, 0.5 * bi)
                        .mul_add(dx, yl)
                        * dx
                };
            } else {
                return 0.0;
            }
        }
        result
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_cspline() {
        let data_x = vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
        let data_y = vec![
            1.0,
            0.961538461538461,
            0.862068965517241,
            0.735294117647059,
            0.609756097560976,
            0.500000000000000,
        ];
        let test_x = vec![
            0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26,
            0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54,
            0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82,
            0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        ];
        let test_y = vec![
            1.000000000000000,
            0.997583282975581,
            0.995079933416512,
            0.992403318788142,
            0.989466806555819,
            0.986183764184894,
            0.982467559140716,
            0.978231558888635,
            0.973389130893999,
            0.967853642622158,
            0.961538461538461,
            0.954382579685350,
            0.946427487413627,
            0.937740299651188,
            0.928388131325928,
            0.918438097365742,
            0.907957312698524,
            0.897012892252170,
            0.885671950954575,
            0.874001603733634,
            0.862068965517241,
            0.849933363488199,
            0.837622973848936,
            0.825158185056786,
            0.812559385569085,
            0.799846963843167,
            0.787041308336369,
            0.774162807506023,
            0.761231849809467,
            0.748268823704033,
            0.735294117647059,
            0.722328486073082,
            0.709394147325463,
            0.696513685724764,
            0.683709685591549,
            0.671004731246381,
            0.658421407009825,
            0.645982297202442,
            0.633709986144797,
            0.621627058157454,
            0.609756097560976,
            0.598112015427308,
            0.586679029833925,
            0.575433685609685,
            0.564352527583445,
            0.553412100584061,
            0.542588949440392,
            0.531859618981294,
            0.521200654035625,
            0.510588599432241,
        ];
        let test_dy = vec![
            -0.120113913432180,
            -0.122279726798445,
            -0.128777166897241,
            -0.139606233728568,
            -0.154766927292426,
            -0.174259247588814,
            -0.198083194617734,
            -0.226238768379184,
            -0.258725968873165,
            -0.295544796099676,
            -0.336695250058719,
            -0.378333644186652,
            -0.416616291919835,
            -0.451543193258270,
            -0.483114348201955,
            -0.511329756750890,
            -0.536189418905076,
            -0.557693334664512,
            -0.575841504029200,
            -0.590633926999137,
            -0.602070603574326,
            -0.611319695518765,
            -0.619549364596455,
            -0.626759610807396,
            -0.632950434151589,
            -0.638121834629033,
            -0.642273812239728,
            -0.645406366983674,
            -0.647519498860871,
            -0.648613207871319,
            -0.648687494015019,
            -0.647687460711257,
            -0.645558211379322,
            -0.642299746019212,
            -0.637912064630930,
            -0.632395167214473,
            -0.625749053769843,
            -0.617973724297039,
            -0.609069178796061,
            -0.599035417266910,
            -0.587872439709585,
            -0.576731233416743,
            -0.566762785681043,
            -0.557967096502484,
            -0.550344165881066,
            -0.543893993816790,
            -0.538616580309654,
            -0.534511925359660,
            -0.531580028966807,
            -0.529820891131095,
        ];
        let test_iy = vec![
            0.000000000000000,
            0.019975905023535,
            0.039902753768792,
            0.059777947259733,
            0.079597153869625,
            0.099354309321042,
            0.119041616685866,
            0.138649546385285,
            0.158166836189794,
            0.177580491219196,
            0.196875783942601,
            0.216036382301310,
            0.235045759060558,
            0.253888601161251,
            0.272550937842853,
            0.291020140643388,
            0.309284923399436,
            0.327335342246135,
            0.345162795617181,
            0.362760024244829,
            0.380121111159890,
            0.397241442753010,
            0.414117280448683,
            0.430745332379281,
            0.447122714446318,
            0.463246950320456,
            0.479115971441505,
            0.494728117018421,
            0.510082134029305,
            0.525177177221407,
            0.540012809111123,
            0.554589001813881,
            0.568906157172889,
            0.582965126887879,
            0.596767214344995,
            0.610314174616794,
            0.623608214462242,
            0.636651992326715,
            0.649448618342004,
            0.662001654326309,
            0.674315113784241,
            0.686393423540581,
            0.698241001711602,
            0.709861835676399,
            0.721259443710643,
            0.732436874986582,
            0.743396709573044,
            0.754141058435429,
            0.764671563435718,
            0.774989397332469,
        ];

        let mut spline = CubicSpline::new(&data_x, &data_y);

        for (i, x) in test_x.iter().enumerate() {
            let y = spline.eval(*x);
            let dy = spline.deriv(*x);
            let iy = spline.integrate(test_x[0], *x);

            let diff_y = y - test_y[i];
            let diff_dy = dy - test_dy[i];
            let diff_iy = iy - test_iy[i];
            println!("{:?}, {:?}", y, test_y[i]);
            println!("{:?}, {:?}", dy, test_dy[i]);
            println!("{:?}, {:?}", iy, test_iy[i]);
            assert!(diff_y.abs() < 1e-10);
            assert!(diff_dy.abs() < 1e-10);
            assert!(diff_iy.abs() < 1e-10);
        }
    }
}
