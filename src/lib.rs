#![forbid(unsafe_code)]

pub(crate) mod details;
pub mod distance;
pub mod fuzz;

#[derive(Debug, Copy, Clone)]
pub enum Hash {
    UNSIGNED(u64),
    SIGNED(i64),
}

pub trait HashableChar {
    fn hash_char(&self) -> Hash;
}
