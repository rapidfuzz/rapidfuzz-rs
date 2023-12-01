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

use crate::common::{DistanceCutoff, NoScoreCutoff, SimilarityCutoff, WithScoreCutoff};
use crate::details::distance::MetricUsize;
use crate::HashableChar;

use std::error;
use std::fmt::{self, Debug, Display, Formatter};

#[derive(Default, Copy, Clone)]
pub struct Padding(bool);
#[derive(Default, Copy, Clone)]
pub struct NoPadding;

#[must_use]
#[derive(Copy, Clone, Debug)]
pub struct Args<ResultType, CutoffType, PaddingType> {
    score_cutoff: CutoffType,
    score_hint: Option<ResultType>,
    pad: PaddingType,
}

impl<ResultType> Default for Args<ResultType, NoScoreCutoff, NoPadding> {
    fn default() -> Args<ResultType, NoScoreCutoff, NoPadding> {
        Args {
            score_cutoff: NoScoreCutoff,
            score_hint: None,
            pad: NoPadding,
        }
    }
}

pub trait PaddingTrait<T>
where
    T: Copy,
{
    type Output: Copy + PartialEq + Debug;

    fn pad(&self) -> bool;
    fn error(&self) -> Self::Output;
    fn score(&self, raw: T) -> Self::Output;
}

impl<T> PaddingTrait<T> for NoPadding
where
    T: Copy + PartialEq + Debug,
{
    type Output = Result<T, Error>;

    fn pad(&self) -> bool {
        false
    }

    fn error(&self) -> Self::Output {
        Err(Error::DifferentLengthArgs)
    }

    fn score(&self, raw: T) -> Self::Output {
        Ok(raw)
    }
}

impl<T> PaddingTrait<T> for Padding
where
    T: Copy + PartialOrd + Debug + Default,
{
    type Output = T;

    fn pad(&self) -> bool {
        self.0
    }

    // will not occur
    fn error(&self) -> Self::Output {
        T::default()
    }

    fn score(&self, raw: T) -> Self::Output {
        raw
    }
}

impl<ResultType, CutoffType, PaddingType> Args<ResultType, CutoffType, PaddingType>
where
    ResultType: Copy,
{
    pub fn score_hint(mut self, score_hint: ResultType) -> Self {
        self.score_hint = Some(score_hint);
        self
    }

    pub fn score_cutoff(
        self,
        score_cutoff: ResultType,
    ) -> Args<ResultType, WithScoreCutoff<ResultType>, PaddingType> {
        Args {
            score_hint: self.score_hint,
            score_cutoff: WithScoreCutoff(score_cutoff),
            pad: self.pad,
        }
    }

    pub fn pad(self, pad: bool) -> Args<ResultType, CutoffType, Padding> {
        Args {
            score_hint: self.score_hint,
            score_cutoff: self.score_cutoff,
            pad: Padding(pad),
        }
    }
}

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

fn distance_impl<Iter1, Iter2>(mut s1: Iter1, mut s2: Iter2) -> usize
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
                return dist;
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
        _score_cutoff: Option<usize>,
        _score_hint: Option<usize>,
    ) -> usize
    where
        Iter1: Iterator,
        Iter2: Iterator,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar,
    {
        distance_impl(s1, s2)
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
/// assert_eq!(Ok(1), hamming::distance("hamming".chars(), "humming".chars()));
/// ```
pub fn distance<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> Result<usize, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    distance_with_args(s1, s2, &Args::default())
}

pub fn distance_with_args<Iter1, Iter2, CutoffType, PaddingType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<usize, CutoffType, PaddingType>,
) -> PaddingType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: DistanceCutoff<usize>,
    PaddingType: PaddingTrait<CutoffType::Output>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !args.pad.pad() && len1 != len2 {
        return args.pad.error();
    }

    args.pad
        .score(args.score_cutoff.score(IndividualComparator {}._distance(
            s1_iter,
            len1,
            s2_iter,
            len2,
            args.score_cutoff.cutoff(),
            args.score_hint,
        )))
}

/// Hamming similarity in the range [0, max]
///
/// This is calculated as `max(len1, len2) - `[`distance`].
///
pub fn similarity<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> Result<usize, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    similarity_with_args(s1, s2, &Args::default())
}

pub fn similarity_with_args<Iter1, Iter2, CutoffType, PaddingType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<usize, CutoffType, PaddingType>,
) -> PaddingType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: SimilarityCutoff<usize>,
    PaddingType: PaddingTrait<CutoffType::Output>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !args.pad.pad() && len1 != len2 {
        return args.pad.error();
    }

    args.pad
        .score(args.score_cutoff.score(IndividualComparator {}._similarity(
            s1_iter,
            len1,
            s2_iter,
            len2,
            args.score_cutoff.cutoff(),
            args.score_hint,
        )))
}

/// Normalized Hamming distance in the range [1.0, 0.0]
///
/// This is calculated as [`distance`]` / max(len1, len2)`.
///
pub fn normalized_distance<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> Result<f64, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    normalized_distance_with_args(s1, s2, &Args::default())
}

pub fn normalized_distance_with_args<Iter1, Iter2, CutoffType, PaddingType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType, PaddingType>,
) -> PaddingType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: DistanceCutoff<f64>,
    PaddingType: PaddingTrait<CutoffType::Output>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !args.pad.pad() && len1 != len2 {
        return args.pad.error();
    }

    args.pad.score(
        args.score_cutoff
            .score(IndividualComparator {}._normalized_distance(
                s1_iter,
                len1,
                s2_iter,
                len2,
                args.score_cutoff.cutoff(),
                args.score_hint,
            )),
    )
}

/// Normalized Hamming similarity in the range [0.0, 1.0]
///
/// This is calculated as `1.0 - `[`normalized_distance`].
///
pub fn normalized_similarity<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> Result<f64, Error>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    normalized_similarity_with_args(s1, s2, &Args::default())
}

pub fn normalized_similarity_with_args<Iter1, Iter2, CutoffType, PaddingType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType, PaddingType>,
) -> PaddingType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: SimilarityCutoff<f64>,
    PaddingType: PaddingTrait<CutoffType::Output>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    if !args.pad.pad() && len1 != len2 {
        return args.pad.error();
    }

    args.pad.score(
        args.score_cutoff
            .score(IndividualComparator {}._normalized_similarity(
                s1_iter,
                len1,
                s2_iter,
                len2,
                args.score_cutoff.cutoff(),
                args.score_hint,
            )),
    )
}

/// `One x Many` comparisons using the Hamming distance
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::hamming;
///
/// let scorer = hamming::BatchComparator::new("hamming".chars());
/// assert_eq!(Ok(1), scorer.distance("humming".chars()));
/// ```
pub struct BatchComparator<Elem1> {
    s1: Vec<Elem1>,
}

impl<Elem1> BatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
    {
        Self {
            s1: s1.into_iter().collect(),
        }
    }

    /// Distance calculated similar to [`distance`]
    pub fn distance<Iter2>(&self, s2: Iter2) -> Result<usize, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        distance(self.s1.iter().copied(), s2)
    }

    pub fn distance_with_args<Iter2, CutoffType, PaddingType>(
        &self,
        s2: Iter2,
        args: &Args<usize, CutoffType, PaddingType>,
    ) -> PaddingType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: DistanceCutoff<usize>,
        PaddingType: PaddingTrait<CutoffType::Output>,
    {
        distance_with_args(self.s1.iter().copied(), s2, args)
    }

    /// Similarity calculated similar to [`similarity`]
    pub fn similarity<Iter2>(&self, s2: Iter2) -> Result<usize, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        similarity(self.s1.iter().copied(), s2)
    }

    pub fn similarity_with_args<Iter2, CutoffType, PaddingType>(
        &self,
        s2: Iter2,
        args: &Args<usize, CutoffType, PaddingType>,
    ) -> PaddingType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: SimilarityCutoff<usize>,
        PaddingType: PaddingTrait<CutoffType::Output>,
    {
        similarity_with_args(self.s1.iter().copied(), s2, args)
    }

    /// Normalized distance calculated similar to [`normalized_distance`]
    pub fn normalized_distance<Iter2>(&self, s2: Iter2) -> Result<f64, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        normalized_distance(self.s1.iter().copied(), s2)
    }

    pub fn normalized_distance_with_args<Iter2, CutoffType, PaddingType>(
        &self,
        s2: Iter2,
        args: &Args<f64, CutoffType, PaddingType>,
    ) -> PaddingType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: DistanceCutoff<f64>,
        PaddingType: PaddingTrait<CutoffType::Output>,
    {
        normalized_distance_with_args(self.s1.iter().copied(), s2, args)
    }

    /// Normalized similarity calculated similar to [`normalized_similarity`]
    pub fn normalized_similarity<Iter2>(&self, s2: Iter2) -> Result<f64, Error>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        normalized_similarity(self.s1.iter().copied(), s2)
    }

    pub fn normalized_similarity_with_args<Iter2, CutoffType, PaddingType>(
        &self,
        s2: Iter2,
        args: &Args<f64, CutoffType, PaddingType>,
    ) -> PaddingType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: SimilarityCutoff<f64>,
        PaddingType: PaddingTrait<CutoffType::Output>,
    {
        normalized_similarity_with_args(self.s1.iter().copied(), s2, args)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_dist(dist: usize, str1: &str, str2: &str) {
        assert_eq!(Ok(dist), distance(str1.chars(), str2.chars()));
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
        assert_eq!(Ok(1), distance([1, 2, 4], [1, 2, 3]));
    }

    #[test]
    fn diff() {
        assert_dist(3, "hamming", "hammers");

        assert_eq!(
            3,
            distance_with_args(
                "hammers".chars(),
                "hamming".chars(),
                &Args::default().pad(true)
            )
        );
        assert_eq!(
            Some(3),
            distance_with_args(
                "hammers".chars(),
                "hamming".chars(),
                &Args::default().pad(true).score_cutoff(3)
            )
        );
        assert_eq!(
            None,
            distance_with_args(
                "hammers".chars(),
                "hamming".chars(),
                &Args::default().pad(true).score_cutoff(2)
            )
        );
        assert_eq!(
            Ok(Some(3)),
            distance_with_args(
                "hammers".chars(),
                "hamming".chars(),
                &Args::default().score_cutoff(3)
            )
        );
        assert_eq!(
            Ok(None),
            distance_with_args(
                "hammers".chars(),
                "hamming".chars(),
                &Args::default().score_cutoff(2)
            )
        );
    }

    #[test]
    fn diff_multibyte() {
        assert_dist(2, "hamming", "h香mmüng");
    }

    #[test]
    fn unequal_length() {
        assert_eq!(
            Err(Error::DifferentLengthArgs),
            distance("ham".chars(), "hamming".chars())
        );

        assert_eq!(
            4,
            distance_with_args("ham".chars(), "hamming".chars(), &Args::default().pad(true))
        );

        assert_eq!(
            None,
            distance_with_args(
                "ham".chars(),
                "hamming".chars(),
                &Args::default().pad(true).score_cutoff(3)
            )
        );
    }

    #[test]
    fn names() {
        assert_dist(14, "Friedrich Nietzs", "Jean-Paul Sartre")
    }
}
