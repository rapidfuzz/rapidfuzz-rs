//! `RapidFuzz` is a general purpose string matching library with implementations
//! for Rust, C++ and Python.
//!
//! ## Key Features
//!
//! - **Diverse String Metrics**: Offers a variety of string metrics
//!   to suit different use cases. These range from the Levenshtein
//!   distance for edit-based comparisions to the Jaro-Winkler similarity for
//!   more nuanced similarity assessments.
//! - **Optimized for Speed**: The library is designed with performance in mind.
//!   Each implementation is carefully designed to ensure optimal performance,
//!   making it suitable for the analysis of large datasets.
//! - **Easy to use**: The API is designed to be simple to use, while still giving
//!   the implementation room for optimization.
//!
//! ## Installation
//!
//! The installation is as simple as:
//! ```console
//! $ cargo add rapidfuzz
//! ```
//!
//! ## Usage
//!
//! The following examples show the usage with the [`Levenshtein`] distance. Other metrics
//! can be found in the [`fuzz`] and [`distance`] modules.
//!
//! ```rust
//! use rapidfuzz::distance::levenshtein;
//!
//! // Perform a simple comparision using he levenshtein distance
//! assert_eq!(
//!     Some(3),
//!     levenshtein::distance("kitten".chars(), "sitting".chars(), None, None, None)
//! );
//!
//! // If you are sure the input strings are ascii only it's usually faster to operate on bytes
//! assert_eq!(
//!     Some(3),
//!     levenshtein::distance("kitten".bytes(), "sitting".bytes(), None, None, None)
//! );
//!
//! // You can provide a score_cutoff value to filter out strings with distance that is worse than
//! // the score_cutoff
//! assert_eq!(
//!     None,
//!     levenshtein::distance("kitten".chars(), "sitting".chars(), None, 2, None)
//! );
//!
//! // You can provide a score_hint to tell the implementation about the expected score.
//! // This can be used to select a more performant implementation internally, but might cause
//! // a slowdown in cases where the distance is actually worse than the score_hint
//! assert_eq!(
//!     Some(3),
//!     levenshtein::distance("kitten".chars(), "sitting".chars(), None, None, Some(3))
//! );
//!
//! // When comparing a single string to multiple strings you can use the provided `BatchComparators`.
//! // These can cache part of the calculation which can provide significant speedups
//! let scorer = levenshtein::BatchComparator::new("kitten".chars(), None);
//! assert_eq!(Some(3), scorer.distance("sitting".chars(), None, None));
//! assert_eq!(Some(0), scorer.distance("kitten".chars(), None, None));
//! ```
//!
//! [`Levenshtein`]: distance/levenshtein/index.html
//! [`fuzz`]: fuzz/index.html
//! [`distance`]: distance/index.html

#![forbid(unsafe_code)]
#![allow(
    // these casts are sometimes needed. They restrict the length of input iterators
    // but there isn't really any way around this except for always working with
    // 128 bit types
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    // things are often more readable this way
    clippy::module_name_repetitions,
    // not practical
    clippy::needless_pass_by_value,
    clippy::similar_names,
    clippy::too_many_lines,
    // noisy
    clippy::missing_errors_doc,
)]

pub(crate) mod details;
pub mod distance;
pub mod fuzz;

/// Hash value in the range `i64::MIN` - `u64::MAX`
#[derive(Debug, Copy, Clone)]
pub enum Hash {
    UNSIGNED(u64),
    SIGNED(i64),
}

/// trait used to map between element types and unique hash values
///
/// `RapidFuzz` already implements this trait for most primitive types.
/// For custom types this trat can be used to support the internal hashmaps.
/// There are a couple of things to keep in mind when implementing this trait:
/// - hashes have to be a unique value in the range `i64::MIN` - `u64::MAX`.
///   If two distinct objects produce the same hash, they will be assumed to be similar
///   by the hashmap.
/// - the hash function should be very fast. For primitive types it can just be the identity
///   function
/// - the hashmaps are optimized for extended ascii, so values in the range 0-255 generally
///   provide a better performance.
///
/// # Example
/// ```
/// use rapidfuzz::distance;
/// use rapidfuzz::{Hash, HashableChar};
///
/// #[derive(PartialEq)]
/// struct MyType {
///     val: u64,
/// }
///
/// impl HashableChar for &MyType {
///     fn hash_char(&self) -> Hash {
///         Hash::UNSIGNED(self.val)
///     }
/// }
///
/// assert_eq!(
///     Some(1),
///     distance::levenshtein::distance(
///         &[MyType { val: 1 }, MyType { val: 1 }],
///         &[MyType { val: 2 }, MyType { val: 1 }],
///         None,
///         None,
///         None,
///     )
/// );
/// ```
pub trait HashableChar {
    fn hash_char(&self) -> Hash;
}
