use crate::{Hash, HashableChar};
use std::iter::{Skip, Take};

pub fn norm_sim_to_norm_dist(score_cutoff: f64) -> f64 {
    let imprecision = 0.00001;
    (1.0 - score_cutoff + imprecision).min(1.0)
}

macro_rules! impl_hashable_char {
    ($base_type:ty, $kind:tt $(, $t:ty)*) => {
        impl HashableChar for $base_type {
            #[inline]
            fn hash_char(&self) -> Hash
            {
                Hash::$kind(*self $(as $t)*)
            }
        }

        impl HashableChar for &$base_type {
            #[inline]
            fn hash_char(&self) -> Hash
            {
                Hash::$kind(**self $(as $t)*)
            }
        }
    }
}

impl_hashable_char!(char, UNSIGNED, u32, u64);
impl_hashable_char!(i8, SIGNED, i64);
impl_hashable_char!(i16, SIGNED, i64);
impl_hashable_char!(i32, SIGNED, i64);
impl_hashable_char!(i64, SIGNED, i64);
impl_hashable_char!(u8, UNSIGNED, u64);
impl_hashable_char!(u16, UNSIGNED, u64);
impl_hashable_char!(u32, UNSIGNED, u64);
impl_hashable_char!(u64, UNSIGNED, u64);

pub fn find_common_prefix<Iter1, Iter2, Elem1, Elem2>(s1: Iter1, s2: Iter2) -> usize
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2>,
    Elem2: PartialEq<Elem1>,
{
    s1.zip(s2)
        .take_while(|(a_char, b_char)| a_char == b_char)
        .count()
}

pub fn find_common_suffix<Iter1, Iter2, Elem1, Elem2>(s1: Iter1, s2: Iter2) -> usize
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2>,
    Elem2: PartialEq<Elem1>,
{
    s1.rev()
        .zip(s2.rev())
        .take_while(|(a_char, b_char)| a_char == b_char)
        .count()
}

pub struct RemovedAffix<Iter1, Iter2, Elem1, Elem2>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2>,
    Elem2: PartialEq<Elem1>,
{
    pub s1: Skip<Take<Iter1>>,
    pub len1: usize,
    pub s2: Skip<Take<Iter2>>,
    pub len2: usize,
    pub prefix_len: usize,
    pub suffix_len: usize,
}

pub fn remove_common_affix<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    mut len1: usize,
    s2: Iter2,
    mut len2: usize,
) -> RemovedAffix<Iter1, Iter2, Elem1, Elem2>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
{
    let suffix_len = find_common_suffix(s1.clone(), s2.clone());
    let s1_iter_no_suffix = s1.take(len1 - suffix_len);
    let s2_iter_no_suffix = s2.take(len2 - suffix_len);
    let prefix_len = find_common_prefix(s1_iter_no_suffix.clone(), s2_iter_no_suffix.clone());
    let s1_iter = s1_iter_no_suffix.skip(prefix_len);
    let s2_iter = s2_iter_no_suffix.skip(prefix_len);
    len1 -= prefix_len + suffix_len;
    len2 -= prefix_len + suffix_len;

    RemovedAffix {
        s1: s1_iter,
        len1,
        s2: s2_iter,
        len2,
        prefix_len,
        suffix_len,
    }
}
