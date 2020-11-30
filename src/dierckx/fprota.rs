/// Apply a Givens rotation to `a` and `b`.
#[inline]
pub(super) fn fprota(cos: f64, sin: f64, a: &mut f64, b: &mut f64) {
    let stor1 = *a;
    let stor2 = *b;
    *b = cos * stor2 + sin * stor1;
    *a = cos * stor1 - sin * stor2;
}
