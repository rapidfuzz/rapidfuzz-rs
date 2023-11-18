use crate::details::common::{norm_sim_to_norm_dist, HashableChar, UnrefIterator};
use crate::details::distance::{
    build_cached_distance_metric_funcs, build_cached_normalized_metric_funcs,
    build_distance_metric_funcs, build_normalized_metric_funcs,
};

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

fn hamming_impl<Iter1, Iter2, Elem1, Elem2>(
    mut s1: Iter1,
    mut s2: Iter2,
    score_cutoff: usize,
) -> usize
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

struct Hamming {}

impl Hamming {
    build_distance_metric_funcs!(Hamming, usize, 0, usize::MAX);

    fn maximum(len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn distance<Iter1, Iter2, Elem1, Elem2>(
        s1: Iter1,
        _len1: usize,
        s2: Iter2,
        _len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1>,
        Iter2: Iterator<Item = Elem2>,
        Elem1: PartialEq<Elem2> + HashableChar,
        Elem2: PartialEq<Elem1> + HashableChar,
    {
        hamming_impl(s1, s2, score_cutoff)
    }
}

pub fn distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<usize, HammingError>
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

    Ok(Hamming::distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into().unwrap_or(usize::MAX),
        score_hint.into().unwrap_or(usize::MAX),
    ))
}

pub fn similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<usize, HammingError>
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

    Ok(Hamming::similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into().unwrap_or(0),
        score_hint.into().unwrap_or(0),
    ))
}

pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<f64, HammingError>
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

    Ok(Hamming::normalized_distance(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into().unwrap_or(1.0),
        score_hint.into().unwrap_or(1.0),
    ))
}

pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    pad: bool,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Result<f64, HammingError>
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

    Ok(Hamming::normalized_similarity(
        s1_iter,
        len1,
        s2_iter,
        len2,
        score_cutoff.into().unwrap_or(0.0),
        score_hint.into().unwrap_or(0.0),
    ))
}

struct CachedHammingImpl<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub(crate) s1: Vec<Elem1>,
}

impl<Elem1> CachedHammingImpl<Elem1>
where
    Elem1: HashableChar + Clone,
{
    build_cached_distance_metric_funcs!(CachedHammingImpl, usize, 0, usize::MAX);

    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();
        CachedHammingImpl { s1 }
    }

    fn maximum(&self, len2: usize) -> usize {
        self.s1.len() + len2
    }

    fn _distance<Iter2, Elem2>(
        &self,
        s2: Iter2,
        _len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter2: Iterator<Item = Elem2>,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        hamming_impl(
            UnrefIterator {
                seq: self.s1.iter(),
            },
            s2,
            score_cutoff,
        )
    }
}

pub struct CachedHamming<Elem1>
where
    Elem1: HashableChar + Clone,
{
    scorer: CachedHammingImpl<Elem1>,
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
        CachedHamming {
            scorer: CachedHammingImpl::new(s1),
        }
    }

    pub fn distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<usize, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        let s2_iter = s2.into_iter();
        let len1 = self.scorer.s1.len();
        let len2 = s2_iter.clone().count();

        if !pad && len1 != len2 {
            return Err(HammingError::DifferentLengthArgs);
        }

        Ok(self.scorer._distance(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(usize::MAX),
            score_hint.into().unwrap_or(usize::MAX),
        ))
    }

    pub fn similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<usize, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        let s2_iter = s2.into_iter();
        let len1 = self.scorer.s1.len();
        let len2 = s2_iter.clone().count();

        if !pad && len1 != len2 {
            return Err(HammingError::DifferentLengthArgs);
        }

        Ok(self.scorer._similarity(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(0),
            score_hint.into().unwrap_or(0),
        ))
    }

    pub fn normalized_distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<f64, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        let len1 = self.scorer.s1.len();
        let len2 = s2_iter.clone().count();

        if !pad && len1 != len2 {
            return Err(HammingError::DifferentLengthArgs);
        }

        Ok(self.scorer._normalized_distance(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(1.0),
            score_hint.into().unwrap_or(1.0),
        ))
    }

    pub fn normalized_similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        pad: bool,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Result<f64, HammingError>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        let len1 = self.scorer.s1.len();
        let len2 = s2_iter.clone().count();

        if !pad && len1 != len2 {
            return Err(HammingError::DifferentLengthArgs);
        }

        Ok(self.scorer._normalized_similarity(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(0.0),
            score_hint.into().unwrap_or(0.0),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_hamming_dist(dist: usize, str1: &str, str2: &str) {
        assert_eq!(
            Ok(dist),
            distance(str1.chars(), str2.chars(), false, None, None)
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
        assert_eq!(Ok(1), distance([1, 2, 4], [1, 2, 3], false, None, None));
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
            distance("ham".chars(), "hamming".chars(), false, None, None)
        );

        assert_eq!(
            Ok(4),
            distance("ham".chars(), "hamming".chars(), true, None, None)
        );
    }

    #[test]
    fn hamming_names() {
        assert_hamming_dist(14, "Friedrich Nietzs", "Jean-Paul Sartre")
    }
}
