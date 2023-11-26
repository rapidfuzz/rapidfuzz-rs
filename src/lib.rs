#![forbid(unsafe_code)]

pub(crate) mod details;
pub mod distance;
pub mod fuzz;

pub use crate::details::common::{Hash, HashableChar};
