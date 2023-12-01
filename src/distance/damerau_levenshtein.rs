//! Damerau-Levenshtein distance
//!
//! The Damerau-Levenshtein distance measures the minimum number of operations required to
//! transform one string into another, considering four types of elementary edits:
//! `insertions`, `deletions`, `substitutions`, and `transpositions of adjacent characters`.
//! A transposition involves swapping two adjacent characters.
//! It does respect triangle inequality, and is thus a metric distance.
//!
//! It's often used in applications where transpositions are common. An example for this would
//! be typing errors involving adjacent characters.
//!
//! # Differences from Levenshtein distance
//!
//! The Damerau-Levenshtein distance includes transpositions as an additional operation,
//! which allows for the swapping of adjacent characters.
//! The [`Levenshtein`] distance considers insertions, deletions, and substitutions only.
//! The Levenshtein distance is computationally less intensive than the Damerau-Levenshtein distance algorithm,
//! as it involves a simpler set of edit operations.
//!
//! # Differences from Optimal String Alignment (OSA) distance
//!
//! While both the Damerau-Levenshtein and [`OSA`] distance include transpositions,
//! they differ in the treatment of transpositions. OSA treats any transposition as a
//! single operation, regardless of whether the transposed characters are adjacent or not.
//! In contrast, the Damerau-Levenshtein distance specifically allows transpositions of adjacent
//! characters.
//!
//! An example where this leads to different results are the strings `CA` and `ÀBC`
//!
//! ```
//! use rapidfuzz::distance::damerau_levenshtein;
//! use rapidfuzz::distance::osa;
//!
//! assert_eq!(2, damerau_levenshtein::distance("CA".chars(), "ABC".chars()));
//! assert_eq!(3, osa::distance("CA".chars(), "ABC".chars()));
//! ```
//!
//! The handling of transpositions in the OSA distance is simpler, which makes it computationally less intensive.
//!
//! [`Levenshtein`]: ../levenshtein/index.html
//! [`OSA`]: ../osa/index.html
//!
//! # Performance
//!
//! The implementation has a runtime complexity of `O(N*M)` and a memory usage of `O(N+M)`.
//! It's based on the paper
//! `Linear space string correction algorithm using the Damerau-Levenshtein distance`
//! from Chunchun Zhao and Sartaj Sahni
//!
//! ![benchmark results](https://raw.githubusercontent.com/maxbachmann/rapidfuzz-rs/main/rapidfuzz-benches/results/damerau_levenshtein.svg)
//!

use crate::common::{DistanceCutoff, NoScoreCutoff, SimilarityCutoff, WithScoreCutoff};
use crate::details::common::remove_common_affix;
use crate::details::distance::MetricUsize;
use crate::details::growing_hashmap::{GrowingHashmap, HybridGrowingHashmap};
use crate::HashableChar;
use std::cmp::{max, min};
use std::mem;

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

#[derive(Clone, Copy, PartialEq, Eq)]
struct RowId {
    val: isize,
}

impl Default for RowId {
    fn default() -> Self {
        Self { val: -1 }
    }
}

/// based on the paper
/// "Linear space string correction algorithm using the Damerau-Levenshtein distance"
/// from Chunchun Zhao and Sartaj Sahni
///
/// todo in c++ this is templated on an integer type which reduced
/// memory usage depending on string lengths
fn distance_zhao<Iter1, Iter2>(s1: Iter1, len1: usize, s2: Iter2, len2: usize) -> usize
where
    Iter1: Iterator,
    Iter2: Iterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar,
{
    let max_val = max(len1, len2) as isize + 1;

    let mut last_row_id = HybridGrowingHashmap::<RowId> {
        map_unsigned: GrowingHashmap::default(),
        map_signed: GrowingHashmap::default(),
        extended_ascii: [RowId::default(); 256],
    };
    let size = len2 + 2;
    let mut fr = vec![max_val; size];
    let mut r1 = vec![max_val; size];
    let mut r: Vec<isize> = (max_val..=max_val).chain(0..(size - 1) as isize).collect();

    for (i, ch1) in s1.enumerate().map(|(i, ch1)| (i + 1, ch1)) {
        mem::swap(&mut r, &mut r1);
        let mut last_col_id: isize = -1;
        let mut last_i2l1 = r[1];
        r[1] = i as isize;
        let mut t = max_val;

        for (j, ch2) in s2.clone().enumerate().map(|(j, ch2)| (j + 1, ch2)) {
            let diag = r1[j] + isize::from(ch1 != ch2);
            let left = r[j] + 1;
            let up = r1[j + 1] + 1;
            let mut temp = min(diag, min(left, up));

            if ch1 == ch2 {
                last_col_id = j as isize; // last occurence of s1_i
                fr[j + 1] = r1[j - 1]; // save H_k-1,j-2
                t = last_i2l1; // save H_i-2,l-1
            } else {
                let k = last_row_id.get(ch2).val;
                let l = last_col_id;

                if j as isize - l == 1 {
                    let transpose = fr[j + 1] + (i as isize - k);
                    temp = min(temp, transpose);
                } else if i as isize - k == 1 {
                    let transpose = t + (j as isize - l);
                    temp = min(temp, transpose);
                }
            }

            last_i2l1 = r[j + 1];
            r[j + 1] = temp;
        }

        last_row_id.get_mut(ch1).val = i as isize;
    }

    r[len2 + 1] as usize
}

fn distance_impl<Iter1, Iter2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> usize
where
    Iter1: DoubleEndedIterator + Clone,
    Iter2: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar,
{
    if score_cutoff < len1.abs_diff(len2) {
        return usize::MAX;
    }

    let affix = remove_common_affix(s1, len1, s2, len2);
    distance_zhao(affix.s1, affix.len1, affix.s2, affix.len2)
}

struct IndividualComparator;

impl MetricUsize for IndividualComparator {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        max(len1, len2)
    }

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        _score_hint: Option<usize>,
    ) -> usize
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        distance_impl(s1, len1, s2, len2, score_cutoff.unwrap_or(usize::MAX))
    }
}

/// Damerau-Levenshtein distance
///
/// Calculates the Damerau-Levenshtein distance.
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::damerau_levenshtein;
///
/// assert_eq!(2, damerau_levenshtein::distance("CA".chars(), "ABC".chars()));
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

/// Damerau-Levenshtein similarity in the range [0, max]
///
/// This is calculated as `max(len1, len2) - `[`distance`].
///
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

/// Normalized Damerau-Levenshtein distance in the range [1.0, 0.0]
///
/// This is calculated as [`distance`]` / max(len1, len2)`.
///
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

/// Normalized Damerau-Levenshtein similarity in the range [0.0, 1.0]
///
/// This is calculated as `1.0 - `[`normalized_distance`].
///
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

/// `One x Many` comparisons using the Damerau-Levenshtein distance
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::damerau_levenshtein;
///
/// let scorer = damerau_levenshtein::BatchComparator::new("CA".chars());
/// assert_eq!(2, scorer.distance("ABC".chars()));
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
        normalized_distance_with_args(self.s1.iter().copied(), s2, args)
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
        normalized_similarity_with_args(self.s1.iter().copied(), s2, args)
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
        distance_with_args(self.s1.iter().copied(), s2, args)
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
        similarity_with_args(self.s1.iter().copied(), s2, args)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static EMPTY: &str = "";
    static TEST: &str = "aaaa";
    static NO_SUFFIX: &str = "aaa";
    static NO_SUFFIX2: &str = "aaab";
    static SWAPPED1: &str = "abaa";
    static SWAPPED2: &str = "baaa";
    static REPLACE_ALL: &str = "bbbb";

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

    fn _test_distance<Iter1, Iter2, CutoffType>(
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

    fn _test_distance_ascii<CutoffType>(
        s1: &str,
        s2: &str,
        args: &Args<usize, CutoffType>,
    ) -> CutoffType::Output
    where
        CutoffType: DistanceCutoff<usize>,
    {
        let res1 = _test_distance(s1.chars(), s2.chars(), args);
        let res2 = _test_distance(s1.bytes(), s2.bytes(), args);

        assert_eq!(res1, res2);
        res1
    }

    fn _test_normalized_similarity<Iter1, Iter2>(
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

    fn _test_normalized_similarity_ascii(
        s1: &str,
        s2: &str,
        args: &Args<f64, WithScoreCutoff<f64>>,
    ) -> Option<f64> {
        let res1 = _test_normalized_similarity(s1.chars(), s2.chars(), args);
        let res2 = _test_normalized_similarity(s1.bytes(), s2.bytes(), args);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    /// levenshtein calculates empty sequence
    #[test]
    fn damerau_levenshtein_empty() {
        assert_eq!(0, _test_distance_ascii(EMPTY, EMPTY, &Args::default()));
        assert_eq!(4, _test_distance_ascii(TEST, EMPTY, &Args::default()));
    }

    /// levenshtein calculates correct distances
    #[test]
    fn damerau_levenshtein_simple() {
        assert_eq!(0, _test_distance_ascii(TEST, TEST, &Args::default()));
        assert_eq!(1, _test_distance_ascii(TEST, NO_SUFFIX, &Args::default()));
        assert_eq!(1, _test_distance_ascii(TEST, NO_SUFFIX2, &Args::default()));
        assert_eq!(
            1,
            _test_distance_ascii(SWAPPED1, SWAPPED2, &Args::default())
        );
        assert_eq!(4, _test_distance_ascii(TEST, REPLACE_ALL, &Args::default()));
        assert_eq!(2, _test_distance_ascii("CA", "ABC", &Args::default()));

        assert_delta!(
            Some(1.0),
            _test_normalized_similarity_ascii(TEST, TEST, &Args::default().score_cutoff(0.0)),
            0.0001
        );
        assert_delta!(
            Some(0.75),
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX, &Args::default().score_cutoff(0.0)),
            0.0001
        );
        assert_delta!(
            Some(0.75),
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX2, &Args::default().score_cutoff(0.0)),
            0.0001
        );
        assert_delta!(
            Some(0.75),
            _test_normalized_similarity_ascii(
                SWAPPED1,
                SWAPPED2,
                &Args::default().score_cutoff(0.0)
            ),
            0.0001
        );
        assert_delta!(
            Some(0.0),
            _test_normalized_similarity_ascii(
                TEST,
                REPLACE_ALL,
                &Args::default().score_cutoff(0.0)
            ),
            0.0001
        );
    }

    #[test]
    fn unicode() {
        assert_eq!(
            5,
            _test_distance("Иванко".chars(), "Петрунко".chars(), &Args::default())
        );

        assert_eq!(
            10,
            _test_distance("ИвaнкoIvan".chars(), "Петрунко".chars(), &Args::default())
        );
    }
}
