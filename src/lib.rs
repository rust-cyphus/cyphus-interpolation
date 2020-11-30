//! Crate for computing interpolating functions.
//!
//! ## Features
//!
//! ### 1D Interpolators
//! - [`UnivariateSpline`] Port of the fortran code DIERCKX,
//! - [`LinearInterp`] a linear interpolator from the GNU Scientific Library,
//! - [`CubicSpline`] a cubic interpolator from the GNU Scientific Library.

#![allow(dead_code, clippy::too_many_arguments, clippy::many_single_char_names)]

pub mod dierckx;
pub mod interp1d;
pub mod prelude;
pub mod traits;
