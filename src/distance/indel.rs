//! Indel distance
//!
//! The Indel distance is a specialized version of the [`Levenshtein`] distance
//! with only insertions and deletions. It can be calculated from the [`Longest Common Subsequence`].
//!
//! Similar to LCS it's commonly used in Bioinformatics applications like DNA sequence analysis, where insertions
//! and deletions play a crucial role in understanding evolutionary relationships and genetic variations.
//!
//! [`Levenshtein`]: ../levenshtein/index.html
//! [`Longest Common Subsequence`]: ../lcs_seq/index.html
//!
//! ## Performance
//!
//! The implementation has a runtime complexity of `O([K/64]*M)` (with `K = MAX(N, score_cutoff)`) and a memory usage of `O(N)`.
//! It's based on the paper `Bit-Parallel LCS-length Computation Revisited` from Heikki Hyyro
//!
//! ![benchmark results](https://raw.githubusercontent.com/rapidfuzz/rapidfuzz-rs/main/rapidfuzz-benches/results/indel.svg)
//!

use crate::common::{DistanceCutoff, NoScoreCutoff, SimilarityCutoff, WithScoreCutoff};
use crate::details::distance::MetricUsize;
use crate::details::pattern_match_vector::BlockPatternMatchVector;
use crate::distance::lcs_seq;
use crate::HashableChar;

#[must_use]
#[derive(Copy, Clone, Debug)]
pub struct Args<ResultType, CutoffType> {
    pub(crate) score_cutoff: CutoffType,
    pub(crate) score_hint: Option<ResultType>,
}

impl<ResultType> Default for Args<ResultType, NoScoreCutoff> {
    fn default() -> Args<ResultType, NoScoreCutoff> {
        Args {
            score_cutoff: NoScoreCutoff,
            score_hint: None,
        }
    }
}

impl<ResultType, CutoffType> Args<ResultType, CutoffType> {
    pub fn score_hint(mut self, score_hint: ResultType) -> Self {
        self.score_hint = Some(score_hint);
        self
    }

    pub fn score_cutoff(
        self,
        score_cutoff: ResultType,
    ) -> Args<ResultType, WithScoreCutoff<ResultType>> {
        Args {
            score_hint: self.score_hint,
            score_cutoff: WithScoreCutoff(score_cutoff),
        }
    }
}

pub(crate) struct IndividualComparator;

impl MetricUsize for IndividualComparator {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1 + len2
    }

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        score_hint_: Option<usize>,
    ) -> usize
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let score_cutoff = score_cutoff_.unwrap_or(usize::MAX);
        let score_hint = score_hint_.unwrap_or(usize::MAX);

        let maximum = self.maximum(len1, len2);
        let lcs_cutoff = if maximum / 2 >= score_cutoff {
            maximum / 2 - score_cutoff
        } else {
            0
        };
        let lcs_hint = if maximum / 2 >= score_hint {
            maximum / 2 - score_hint
        } else {
            0
        };
        let lcs_sim = lcs_seq::IndividualComparator {}._similarity(
            s1,
            len1,
            s2,
            len2,
            Some(lcs_cutoff),
            Some(lcs_hint),
        );
        maximum - 2 * lcs_sim
    }
}

/// Indel distance
///
/// Calculates the minimum number of insertions and deletions required to change
/// one sequence into the other. This is equivalent to the Levenshtein distance
/// with a substitution weight of 2.
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::indel;
///
/// // Find the Indel distance between two strings
/// assert_eq!(3, indel::distance("lewenstein".chars(), "levenshtein".chars()));
///
/// // Setting a maximum distance allows the implementation to select a more efficient implementation
/// assert_eq!(None, indel::distance_with_args("lewenstein".chars(), "levenshtein".chars(), &indel::Args::default().score_cutoff(2)));
/// ```
pub fn distance<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> usize
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

pub fn distance_with_args<Iter1, Iter2, CutoffType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<usize, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: DistanceCutoff<usize>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    args.score_cutoff.score(IndividualComparator {}._distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        args.score_cutoff.cutoff(),
        args.score_hint,
    ))
}

/// Indel similarity in the range [0, max]
///
/// This is calculated as `(len1 + len2) - `[`distance`].
pub fn similarity<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> usize
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

pub fn similarity_with_args<Iter1, Iter2, CutoffType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<usize, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: SimilarityCutoff<usize>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    args.score_cutoff.score(IndividualComparator {}._similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        args.score_cutoff.cutoff(),
        args.score_hint,
    ))
}

/// Normalized Indel distance in the range [1.0, 0.0]
///
/// This is calculated as [`distance`]` / (len1 + len2)`.
pub fn normalized_distance<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> f64
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

pub fn normalized_distance_with_args<Iter1, Iter2, CutoffType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: DistanceCutoff<f64>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    args.score_cutoff
        .score(IndividualComparator {}._normalized_distance(
            s1_iter.clone(),
            s1_iter.count(),
            s2_iter.clone(),
            s2_iter.count(),
            args.score_cutoff.cutoff(),
            args.score_hint,
        ))
}

/// Normalized Indel similarity in the range [0.0, 1.0]
///
/// This is calculated as `1.0 - `[`normalized_distance`].
pub fn normalized_similarity<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> f64
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

pub fn normalized_similarity_with_args<Iter1, Iter2, CutoffType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: SimilarityCutoff<f64>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    args.score_cutoff
        .score(IndividualComparator {}._normalized_similarity(
            s1_iter.clone(),
            s1_iter.count(),
            s2_iter.clone(),
            s2_iter.count(),
            args.score_cutoff.cutoff(),
            args.score_hint,
        ))
}

pub(crate) fn distance_with_pm<Iter1, Iter2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> usize
where
    Iter1: DoubleEndedIterator + Clone,
    Iter2: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    let maximum = len1 + len2;
    let lcs_cutoff = if maximum / 2 >= score_cutoff {
        maximum / 2 - score_cutoff
    } else {
        0
    };

    let lcs_sim = lcs_seq::similarity_with_pm(pm, s1, len1, s2, len2, lcs_cutoff);
    maximum - 2 * lcs_sim
}

/// `One x Many` comparisons using the Indel distance
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::indel;
///
/// let scorer = indel::BatchComparator::new("lewenstein".chars());
/// assert_eq!(3, scorer.distance("levenshtein".chars()));
/// ```
#[derive(Clone)]
pub struct BatchComparator<Elem1> {
    pub(crate) scorer: lcs_seq::BatchComparator<Elem1>,
}

impl<CharT> MetricUsize for BatchComparator<CharT>
where
    CharT: HashableChar + Clone,
{
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1 + len2
    }

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        score_hint_: Option<usize>,
    ) -> usize
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let score_cutoff = score_cutoff_.unwrap_or(usize::MAX);
        let score_hint = score_hint_.unwrap_or(usize::MAX);

        let maximum = self.maximum(len1, len2);
        let lcs_cutoff = if maximum / 2 >= score_cutoff {
            maximum / 2 - score_cutoff
        } else {
            0
        };
        let lcs_hint = if maximum / 2 >= score_hint {
            maximum / 2 - score_hint
        } else {
            0
        };
        let lcs_sim = self
            .scorer
            ._similarity(s1, len1, s2, len2, Some(lcs_cutoff), Some(lcs_hint));
        maximum - 2 * lcs_sim
    }
}

impl<Elem1> BatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        Self {
            scorer: lcs_seq::BatchComparator::new(s1),
        }
    }

    /// Normalized distance calculated similar to [`normalized_distance`]
    pub fn normalized_distance<Iter2>(&self, s2: Iter2) -> f64
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.normalized_distance_with_args(s2, &Args::default())
    }

    pub fn normalized_distance_with_args<Iter2, CutoffType>(
        &self,
        s2: Iter2,
        args: &Args<f64, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: DistanceCutoff<f64>,
    {
        let s2_iter = s2.into_iter();
        args.score_cutoff.score(self._normalized_distance(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            args.score_cutoff.cutoff(),
            args.score_hint,
        ))
    }

    /// Normalized similarity calculated similar to [`normalized_similarity`]
    pub fn normalized_similarity<Iter2>(&self, s2: Iter2) -> f64
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.normalized_similarity_with_args(s2, &Args::default())
    }

    pub fn normalized_similarity_with_args<Iter2, CutoffType>(
        &self,
        s2: Iter2,
        args: &Args<f64, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: SimilarityCutoff<f64>,
    {
        let s2_iter = s2.into_iter();
        args.score_cutoff.score(self._normalized_similarity(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            args.score_cutoff.cutoff(),
            args.score_hint,
        ))
    }

    /// Distance calculated similar to [`distance`]
    pub fn distance<Iter2>(&self, s2: Iter2) -> usize
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.distance_with_args(s2, &Args::default())
    }

    pub fn distance_with_args<Iter2, CutoffType>(
        &self,
        s2: Iter2,
        args: &Args<usize, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: DistanceCutoff<usize>,
    {
        let s2_iter = s2.into_iter();
        args.score_cutoff.score(self._distance(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            args.score_cutoff.cutoff(),
            args.score_hint,
        ))
    }

    /// Similarity calculated similar to [`similarity`]
    pub fn similarity<Iter2>(&self, s2: Iter2) -> usize
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.similarity_with_args(s2, &Args::default())
    }

    pub fn similarity_with_args<Iter2, CutoffType>(
        &self,
        s2: Iter2,
        args: &Args<usize, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: SimilarityCutoff<usize>,
    {
        let s2_iter = s2.into_iter();
        args.score_cutoff.score(self._similarity(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            args.score_cutoff.cutoff(),
            args.score_hint,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            match ($x, $y) {
                (None, None) => {}
                (Some(val1), Some(val2)) => {
                    if (val1 - val2).abs() > $d {
                        panic!("{:?} != {:?}", $x, $y);
                    }
                }
                (_, _) => panic!("{:?} != {:?}", $x, $y),
            }
        };
    }

    fn test_distance<Iter1, Iter2, CutoffType>(
        s1_: Iter1,
        s2_: Iter2,
        args: &Args<usize, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
        CutoffType: DistanceCutoff<usize>,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = distance_with_args(s1.clone(), s2.clone(), args);
        let res2 = distance_with_args(s2.clone(), s1.clone(), args);

        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.distance_with_args(s2.clone(), args);
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.distance_with_args(s1.clone(), args);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn test_distance_ascii<CutoffType>(
        s1: &str,
        s2: &str,
        args: &Args<usize, CutoffType>,
    ) -> CutoffType::Output
    where
        CutoffType: DistanceCutoff<usize>,
    {
        let res1 = test_distance(s1.chars(), s2.chars(), args);
        let res2 = test_distance(s1.bytes(), s2.bytes(), args);

        assert_eq!(res1, res2);
        res1
    }

    fn test_similarity<Iter1, Iter2, CutoffType>(
        s1_: Iter1,
        s2_: Iter2,
        args: &Args<usize, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
        CutoffType: SimilarityCutoff<usize>,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = similarity_with_args(s1.clone(), s2.clone(), args);
        let res2 = similarity_with_args(s2.clone(), s1.clone(), args);

        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.similarity_with_args(s2.clone(), args);
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.similarity_with_args(s1.clone(), args);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn test_similarity_ascii<CutoffType>(
        s1: &str,
        s2: &str,
        args: &Args<usize, CutoffType>,
    ) -> CutoffType::Output
    where
        CutoffType: SimilarityCutoff<usize>,
    {
        let res1 = test_similarity(s1.chars(), s2.chars(), args);
        let res2 = test_similarity(s1.bytes(), s2.bytes(), args);

        assert_eq!(res1, res2);
        res1
    }

    fn test_normalized_distance<Iter1, Iter2>(
        s1_: Iter1,
        s2_: Iter2,
        args: &Args<f64, WithScoreCutoff<f64>>,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = normalized_distance_with_args(s1.clone(), s2.clone(), args);
        let res2 = normalized_distance_with_args(s2.clone(), s1.clone(), args);
        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.normalized_distance_with_args(s2.clone(), args);
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.normalized_distance_with_args(s1.clone(), args);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn test_normalized_distance_ascii(
        s1: &str,
        s2: &str,
        args: &Args<f64, WithScoreCutoff<f64>>,
    ) -> Option<f64> {
        let res1 = test_normalized_distance(s1.chars(), s2.chars(), args);
        let res2 = test_normalized_distance(s1.bytes(), s2.bytes(), args);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    fn test_normalized_similarity<Iter1, Iter2>(
        s1_: Iter1,
        s2_: Iter2,
        args: &Args<f64, WithScoreCutoff<f64>>,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = normalized_similarity_with_args(s1.clone(), s2.clone(), args);
        let res2 = normalized_similarity_with_args(s2.clone(), s1.clone(), args);
        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.normalized_similarity_with_args(s2.clone(), args);
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.normalized_similarity_with_args(s1.clone(), args);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn test_normalized_similarity_ascii(
        s1: &str,
        s2: &str,
        args: &Args<f64, WithScoreCutoff<f64>>,
    ) -> Option<f64> {
        let res1 = test_normalized_similarity(s1.chars(), s2.chars(), args);
        let res2 = test_normalized_similarity(s1.bytes(), s2.bytes(), args);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    #[test]
    fn similar() {
        assert_eq!(0, test_distance_ascii("aaaa", "aaaa", &Args::default()));
        assert_eq!(8, test_similarity_ascii("aaaa", "aaaa", &Args::default()));
        assert_delta!(
            Some(0.0),
            test_normalized_distance_ascii("aaaa", "aaaa", &Args::default().score_cutoff(1.0)),
            0.0001
        );
        assert_delta!(
            Some(1.0),
            test_normalized_similarity_ascii("aaaa", "aaaa", &Args::default().score_cutoff(0.0)),
            0.0001
        );
    }

    #[test]
    fn completely_different() {
        assert_eq!(8, test_distance_ascii("aaaa", "bbbb", &Args::default()));
        assert_eq!(0, test_similarity_ascii("aaaa", "bbbb", &Args::default()));
        assert_delta!(
            Some(1.0),
            test_normalized_distance_ascii("aaaa", "bbbb", &Args::default().score_cutoff(1.0)),
            0.0001
        );
        assert_delta!(
            Some(0.0),
            test_normalized_similarity_ascii("aaaa", "bbbb", &Args::default().score_cutoff(0.0)),
            0.0001
        );
    }

    #[test]
    fn test_mbleven() {
        let mut a = "South Korea";
        let mut b = "North Korea";

        assert_eq!(4, test_distance_ascii(a, b, &Args::default()));
        assert_eq!(
            Some(4),
            test_distance_ascii(a, b, &Args::default().score_cutoff(5))
        );
        assert_eq!(
            Some(4),
            test_distance_ascii(a, b, &Args::default().score_cutoff(4))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(3))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(2))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(1))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(0))
        );

        a = "aabc";
        b = "cccd";
        assert_eq!(6, test_distance_ascii(a, b, &Args::default()));
        assert_eq!(
            Some(6),
            test_distance_ascii(a, b, &Args::default().score_cutoff(6))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(5))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(4))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(3))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(2))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(1))
        );
        assert_eq!(
            None,
            test_distance_ascii(a, b, &Args::default().score_cutoff(0))
        );
    }

    // this was an issue in the cached lcs implementation of rapidfuzz-cpp
    #[test]
    fn test_issue_unknown() {
        let s1 = "001";
        let s2 = "220";
        assert_delta!(
            Some(0.3333333),
            test_normalized_similarity_ascii(s1, s2, &Args::default().score_cutoff(0.0)),
            0.0001
        );
    }

    #[test]
    fn test_banded_implementation() {
        let mut s1 = "ddccbccc";
        let mut s2 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
            aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacca\
            cccaccaaaaaaaadaaaaaaaaccccaccccccaaaaaaaccccaaacccaccccadddaaaaaaaaaaaaaaaaa\
            aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaccccccccacccaaaaaacccaaaaaacc\
            cacccaaaaaacccdccccccaccccccccccccccccccccccccccccccccccccccccccccccccccccccc\
            ccccccddddddaaaaaaaaaaaaaaaaaaaaaaaaaacacccaaaaaacccddddaaaaaaaaaaaaaaaaaaaaa\
            aaaaaaaaccccaaaaaaaaaaccccccaadddaaaaaaaaaaaaaaaaaaaaaacaaaaaa";
        assert_eq!(508, test_distance_ascii(s1, s2, &Args::default()));
        assert_eq!(
            Some(508),
            test_distance_ascii(s1, s2, &Args::default().score_cutoff(508))
        );
        assert_eq!(
            None,
            test_distance_ascii(s1, s2, &Args::default().score_cutoff(507))
        );
        assert_eq!(
            Some(508),
            test_distance_ascii(s1, s2, &Args::default().score_cutoff(usize::MAX))
        );

        s1 = "bbbdbbmbbbbbbbbbBbfbbbbbbbbbbbbbbbbbbbrbbbbbrbbbbbdbnbbbjbhbbbbbbbbbhbbb\
            bbCbobbbxbbbbbkbbbAbxbbwbbbtbcbbbbebbiblbbbbqbbbbbbpbbbbbbubbbkbbDbbbhbkbC\
            bbgbbrbbbbbbbbbbbkbyvbbsbAbbbbz";
        s2 = "jaaagaaqyaaaanrCfwaaxaeahtaaaCzaaaspaaBkvaaaaqDaacndaaeolwiaaauaaaaaaamA";

        assert_eq!(231, test_distance_ascii(s1, s2, &Args::default()));
    }

    #[test]
    fn unicode() {
        assert_eq!(
            8,
            test_distance("Иванко".chars(), "Петрунко".chars(), &Args::default())
        );
    }

    #[test]
    fn fuzzing_regressions() {
        assert_eq!(
            2,
            test_distance("ab".chars(), "ac".chars(), &Args::default())
        );
    }
}
