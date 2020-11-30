/// Compute a givens rotation angles given piv and ww.
pub(super) fn fpgivs(piv: f64, ww: &mut f64, cos: &mut f64, sin: &mut f64) {
    let store = piv.abs();
    let dd = if store >= *ww {
        store * (1.0 + (*ww / piv).powi(2)).sqrt()
    } else {
        *ww * (1.0 + (piv / *ww).powi(2)).sqrt()
    };
    *cos = *ww / dd;
    *sin = piv / dd;
    *ww = dd;
}
