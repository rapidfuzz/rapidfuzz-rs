use crate::details::common::{norm_sim_to_norm_dist, HashableChar};
use crate::details::distance::{build_distance_metric_funcs, build_normalized_metric_funcs};

use std::error::Error;
use std::fmt::{self, Display, Formatter};

#[derive(Debug, PartialEq)]
pub enum HammingError {
    DifferentLengthArgs,
}

impl Display for HammingError {
    fn fmt(&self, fmt: &mut Formatter) -> Result<(), fmt::Error> {
        let text = match self {
            HammingError::DifferentLengthArgs => "Differing length arguments provided",
        };

        write!(fmt, "{}", text)
    }
}

impl Error for HammingError {}

struct Hamming {}

impl Hamming {
    build_distance_metric_funcs!(Hamming, usize, 0, usize::MAX);

    fn maximum<Iter1, Iter2, Elem1, Elem2>(
        _s1: Iter1,
        len1: usize,
        _s2: Iter2,
        len2: usize,
    ) -> usize
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
        Elem1: PartialEq<Elem2> + HashableChar,
        Elem2: PartialEq<Elem1> + HashableChar,
    {
        len1.max(len2)
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        s1: Iter1,
        _len1: usize,
        s2: Iter2,
        _len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
        Elem1: PartialEq<Elem2> + HashableChar,
        Elem2: PartialEq<Elem1> + HashableChar,
    {
        let (mut it_s1, mut it_s2) = (s1.into_iter(), s2.into_iter());
        let mut dist = 0;
        loop {
            match (it_s1.next(), it_s2.next()) {
                (Some(ch1), Some(ch2)) => {
                    if ch1 != ch2 {
                        dist += 1
                    }
                }
                (None, None) => {
                    if dist <= score_cutoff {
                        return dist;
                    } else {
                        return score_cutoff + 1;
                    }
                }
                _ => {
                    dist += 1;
                }
            }
        }
    }
}

pub fn hamming_distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: Option<usize>,
    score_hint: Option<usize>,
) -> Result<usize, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    let score_cutoff = score_cutoff.unwrap_or(usize::MAX);
    let score_hint = score_hint.unwrap_or(usize::MAX);
    Ok(Hamming::_distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff,
        score_hint,
    ))
}

pub fn hamming_similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: Option<usize>,
    score_hint: Option<usize>,
) -> Result<usize, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    let score_cutoff = score_cutoff.unwrap_or(0);
    let score_hint = score_hint.unwrap_or(0);
    Ok(Hamming::_similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff,
        score_hint,
    ))
}

pub fn hamming_normalized_distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: Option<f64>,
    score_hint: Option<f64>,
) -> Result<f64, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    let score_cutoff = score_cutoff.unwrap_or(1.0);
    let score_hint = score_hint.unwrap_or(1.0);
    Ok(Hamming::_normalized_distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff,
        score_hint,
    ))
}

pub fn hamming_normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: Option<f64>,
    score_hint: Option<f64>,
) -> Result<f64, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    let score_cutoff = score_cutoff.unwrap_or(0.0);
    let score_hint = score_hint.unwrap_or(0.0);
    Ok(Hamming::_normalized_similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff,
        score_hint,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_hamming_dist(dist: usize, str1: &str, str2: &str) {
        assert_eq!(
            Ok(dist),
            hamming_distance(str1.chars(), str2.chars(), false, None, None)
        );
    }

    #[test]
    fn hamming_empty() {
        assert_hamming_dist(0, "", "")
    }

    #[test]
    fn hamming_same() {
        assert_hamming_dist(0, "hamming", "hamming")
    }

    #[test]
    fn hamming_numbers() {
        assert_eq!(
            Ok(1),
            hamming_distance([1, 2, 4], [1, 2, 3], false, None, None)
        );
    }

    #[test]
    fn hamming_diff() {
        assert_hamming_dist(3, "hamming", "hammers")
    }

    #[test]
    fn hamming_diff_multibyte() {
        assert_hamming_dist(2, "hamming", "h香mmüng");
    }

    #[test]
    fn hamming_unequal_length() {
        assert_eq!(
            Err(HammingError::DifferentLengthArgs),
            hamming_distance("ham".chars(), "hamming".chars(), false, None, None)
        );

        assert_eq!(
            Ok(4),
            hamming_distance("ham".chars(), "hamming".chars(), true, None, None)
        );
    }

    #[test]
    fn hamming_names() {
        assert_hamming_dist(14, "Friedrich Nietzs", "Jean-Paul Sartre")
    }
}
