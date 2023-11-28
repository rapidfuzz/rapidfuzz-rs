//! Hamming distance
//!
//! The Hamming distance measures the similarity of two sequences of equal length.
//! Specifically, it counts the minimum number of substitutions required to
//! transform one string into the other.
//!
//! While regularly the Hamming distance only works with texts of equal length,
//! this implementation provides an addition argument `pad` to decide whether texts
//! of unequal length should be padded or return an error.
//!

use crate::details::distance::MetricUsize;
use crate::HashableChar;

use std::error;
use std::fmt::{self, Display, Formatter};

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Error {
    DifferentLengthArgs,
}

impl Display for Error {
    fn fmt(&self, fmt: &mut Formatter) -> Result<(), fmt::Error> {
        let text = match self {
            Self::DifferentLengthArgs => "Differing length arguments provided",
        };

        write!(fmt, "{text}")
    }
}

impl error::Error for Error {}

fn distance_impl<Iter1, Iter2>(mut s1: Iter1, mut s2: Iter2, score_cutoff: usize) -> Option<usize>
where
    Iter1: Iterator,
    Iter2: Iterator,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar,
{
    let mut dist = 0;
    loop {
        match (s1.next(), s2.next()) {
            (Some(ch1), Some(ch2)) => {
                if ch1 != ch2 {
                    dist += 1;
                }
            }
            (None, None) => {
                if dist <= score_cutoff {
                    return Some(dist);
                }
                return None;
            }
            _ => {
                dist += 1;
            }
        }
    }
}

struct IndividualComparator;

impl MetricUsize for IndividualComparator {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        _len1: usize,
        s2: Iter2,
        _len2: usize,
        score_cutoff: Option<usize>,
        _score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator,
        Iter2: Iterator,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar,
    {
        distance_impl(s1, s2, score_cutoff.unwrap_or(usize::MAX))
    }
}

/// Hamming distance
///
/// Calculates the Hamming distance.
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::hamming;
///
/// assert_eq!(Ok(Some(1)), hamming::distance("hamming".chars(), "humming".chars(), false, None, None));
/// ```
pub fn distance<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<usize>, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    ScoreCutoff: Into<Option<usize>>,
    ScoreHint: Into<Option<usize>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(Error::DifferentLengthArgs);
    }

    Ok(IndividualComparator {}._distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

/// Hamming similarity in the range [0, max]
///
/// This is calculated as `max(len1, len2) - `[`distance`].
///
pub fn similarity<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<usize>, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    ScoreCutoff: Into<Option<usize>>,
    ScoreHint: Into<Option<usize>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(Error::DifferentLengthArgs);
    }

    Ok(IndividualComparator {}._similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

/// Normalized Hamming distance in the range [1.0, 0.0]
///
/// This is calculated as [`distance`]` / max(len1, len2)`.
///
pub fn normalized_distance<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<f64>, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(Error::DifferentLengthArgs);
    }

    Ok(IndividualComparator {}._normalized_distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

/// Normalized Hamming similarity in the range [0.0, 1.0]
///
/// This is calculated as `1.0 - `[`normalized_distance`].
///
pub fn normalized_similarity<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<f64>, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(Error::DifferentLengthArgs);
    }

    Ok(IndividualComparator {}._normalized_similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

/// `One x Many` comparisons using the Hamming distance
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::hamming;
///
/// let scorer = hamming::BatchComparator::new("hamming".chars(), false);
/// assert_eq!(Ok(Some(1)), scorer.distance("humming".chars(), None, None));
/// ```
pub struct BatchComparator<Elem1> {
    s1: Vec<Elem1>,
    pad: bool,
}

impl<Elem1> BatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1, pad: bool) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
    {
        Self {
            s1: s1.into_iter().collect(),
            pad,
        }
    }

    /// Distance calculated similar to [`distance`]
    pub fn distance<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<usize>, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        distance(
            self.s1.iter().copied(),
            s2,
            self.pad,
            score_cutoff,
            score_hint,
        )
    }

    /// Similarity calculated similar to [`similarity`]
    pub fn similarity<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<usize>, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        similarity(
            self.s1.iter().copied(),
            s2,
            self.pad,
            score_cutoff,
            score_hint,
        )
    }

    /// Normalized distance calculated similar to [`normalized_distance`]
    pub fn normalized_distance<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<f64>, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        normalized_distance(
            self.s1.iter().copied(),
            s2,
            self.pad,
            score_cutoff,
            score_hint,
        )
    }

    /// Normalized similarity calculated similar to [`normalized_similarity`]
    pub fn normalized_similarity<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<f64>, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        normalized_similarity(
            self.s1.iter().copied(),
            s2,
            self.pad,
            score_cutoff,
            score_hint,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_dist(dist: usize, str1: &str, str2: &str) {
        assert_eq!(
            Ok(Some(dist)),
            distance(str1.chars(), str2.chars(), false, None, None)
        );
    }

    #[test]
    fn empty() {
        assert_dist(0, "", "")
    }

    #[test]
    fn same() {
        assert_dist(0, "hamming", "hamming")
    }

    #[test]
    fn numbers() {
        assert_eq!(
            Ok(Some(1)),
            distance([1, 2, 4], [1, 2, 3], false, None, None)
        );
    }

    #[test]
    fn diff() {
        assert_dist(3, "hamming", "hammers")
    }

    #[test]
    fn diff_multibyte() {
        assert_dist(2, "hamming", "h香mmüng");
    }

    #[test]
    fn unequal_length() {
        assert_eq!(
            Err(Error::DifferentLengthArgs),
            distance("ham".chars(), "hamming".chars(), false, None, None)
        );

        assert_eq!(
            Ok(Some(4)),
            distance("ham".chars(), "hamming".chars(), true, None, None)
        );
    }

    #[test]
    fn names() {
        assert_dist(14, "Friedrich Nietzs", "Jean-Paul Sartre")
    }
}
