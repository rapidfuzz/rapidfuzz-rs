#![forbid(unsafe_code)]

pub(crate) mod details;
pub mod distance;
pub mod fuzz;

/// Hash value in the range i64::MIN - u64::MAX
#[derive(Debug, Copy, Clone)]
pub enum Hash {
    UNSIGNED(u64),
    SIGNED(i64),
}

/// trait used to map between element types and unique hash values
///
/// RapidFuzz already implements this trait for most primitive types.
/// For custom types this trat can be used to support the internal hashmaps.
/// There are a couple of things to keep in mind when implementing this trait:
/// - hashes have to be a unique value in the range i64::MIN - u64::MAX.
///   If two distinct objects produce the same hash, they will be assumed to be similar
///   by the hashmap.
/// - the hash function should be very fast. For primitive types it can just be the identity
///   function
/// - the hashmaps are optimized for extended ascii, so values in the range 0-255 generally
///   provide a better performance.
///
pub trait HashableChar {
    fn hash_char(&self) -> Hash;
}
