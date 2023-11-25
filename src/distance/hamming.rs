use crate::details::common::HashableChar;
use crate::details::distance::MetricUsize;

use std::error::Error;
use std::fmt::{self, Display, Formatter};

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum HammingError {
    DifferentLengthArgs,
}

impl Display for HammingError {
    fn fmt(&self, fmt: &mut Formatter) -> Result<(), fmt::Error> {
        let text = match self {
            Self::DifferentLengthArgs => "Differing length arguments provided",
        };

        write!(fmt, "{text}")
    }
}

impl Error for HammingError {}

fn distance_impl<Iter1, Iter2, Elem1, Elem2>(
    mut s1: Iter1,
    mut s2: Iter2,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1>,
    Iter2: Iterator<Item = Elem2>,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
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

struct Hamming;

impl MetricUsize for Hamming {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        _len1: usize,
        s2: Iter2,
        _len2: usize,
        score_cutoff: Option<usize>,
        _score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1>,
        Iter2: Iterator<Item = Elem2>,
        Elem1: PartialEq<Elem2> + HashableChar,
        Elem2: PartialEq<Elem1> + HashableChar,
    {
        distance_impl(s1, s2, score_cutoff.unwrap_or(usize::MAX))
    }
}

pub fn distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<usize>, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    ScoreCutoff: Into<Option<usize>>,
    ScoreHint: Into<Option<usize>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    Ok(Hamming {}._distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

pub fn similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<usize>, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    ScoreCutoff: Into<Option<usize>>,
    ScoreHint: Into<Option<usize>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    Ok(Hamming {}._similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<f64>, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    Ok(Hamming {}._normalized_distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

// todo this api is a bit of an outlier, since it is the only one returning an
// error. Should the return just be Result and store score < score_cutoff as error as well?
// having to unwrap both seems like a pain
pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<Option<f64>, HammingError>
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !pad && len1 != len2 {
        return Err(HammingError::DifferentLengthArgs);
    }

    Ok(Hamming {}._normalized_similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into(),
        score_hint.into(),
    ))
}

pub struct CachedHamming<Elem1>
where
    Elem1: HashableChar + Clone,
{
    s1: Vec<Elem1>,
}

impl<Elem1> CachedHamming<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        Self {
            s1: s1.into_iter().collect(),
        }
    }

    pub fn distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<usize>, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        distance(self.s1.iter().copied(), s2, pad, score_cutoff, score_hint)
    }

    pub fn similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<usize>, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        similarity(self.s1.iter().copied(), s2, pad, score_cutoff, score_hint)
    }

    pub fn normalized_distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<f64>, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        normalized_distance(self.s1.iter().copied(), s2, pad, score_cutoff, score_hint)
    }

    pub fn normalized_similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<Option<f64>, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        normalized_similarity(self.s1.iter().copied(), s2, pad, score_cutoff, score_hint)
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
            Err(HammingError::DifferentLengthArgs),
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
