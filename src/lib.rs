//! Crate for computing interpolating functions.
//!
//! ## Features
//!
//! ### 1D Interpolators
//! - [`UnivariateSpline`] Port of the fortran code DIERCKX,
//! - [`LinearInterp`] a linear interpolator from the GNU Scientific Library,
//! - [`CubicSpline`] a cubic interpolator from the GNU Scientific Library.

pub mod dierckx;
pub mod interp1d;
pub mod prelude;
