pub(crate) mod damerau_levenshtein;
pub(crate) mod hamming;
pub(crate) mod indel;
pub(crate) mod lcs_seq;
pub(crate) mod levenshtein;
pub(crate) mod postfix;
pub(crate) mod prefix;

pub use damerau_levenshtein::{
    damerau_levenshtein_distance, damerau_levenshtein_normalized_distance,
    damerau_levenshtein_normalized_similarity, damerau_levenshtein_similarity,
    CachedDamerauLevenshtein,
};

pub use levenshtein::{
    levenshtein_distance, levenshtein_normalized_distance, levenshtein_normalized_similarity,
    levenshtein_similarity, CachedLevenshtein, LevenshteinWeightTable,
};

pub use indel::{
    indel_distance, indel_normalized_distance, indel_normalized_similarity, indel_similarity,
    CachedIndel,
};

pub use lcs_seq::{
    lcs_seq_distance, lcs_seq_normalized_distance, lcs_seq_normalized_similarity,
    lcs_seq_similarity, CachedLcsSeq,
};

pub use hamming::{
    hamming_distance, hamming_normalized_distance, hamming_normalized_similarity,
    hamming_similarity, HammingError,
};

pub use prefix::{
    prefix_distance, prefix_normalized_distance, prefix_normalized_similarity, prefix_similarity,
    CachedPrefix,
};

pub use postfix::{
    postfix_distance, postfix_normalized_distance, postfix_normalized_similarity,
    postfix_similarity, CachedPostfix,
};

#[cfg(test)]
pub(crate) mod example;
