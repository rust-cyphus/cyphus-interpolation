#[allow(dead_code)]
pub(super) fn fprati(
    p1: &mut f64,
    f1: &mut f64,
    p2: f64,
    f2: f64,
    p3: &mut f64,
    f3: &mut f64,
) -> f64 {
    let p = if *p3 > 0.0 {
        //  value of p in case p3 ^= infinity.
        let h1 = *f1 * (f2 - *f3);
        let h2 = f2 * (*f3 - *f1);
        let h3 = *f3 * (*f1 - f2);
        -(*p1 * p2 * h3 + p2 * *p3 * h1 + *p3 * *p1 * h2) / (*p1 * h1 + p2 * h2 + *p3 * h3)
    } else {
        //  value of p in case p3 = infinity.
        (*p1 * (*f1 - *f3) * f2 - p2 * (f2 - *f3) * *f1) / ((*f1 - f2) * *f3)
    };
    // adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
    if f2 >= 0.0 {
        //goto _30;
        *p1 = p2;
        *f1 = f2;
        return p;
    }
    *p3 = p2;
    *f3 = f2;
    return p;
}
