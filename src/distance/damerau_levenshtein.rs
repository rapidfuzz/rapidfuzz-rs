//! Damerau-Levenshtein distance
//!
//! Similar to [`Levenshtein`], the Damerau-Levenshtein distance is the minimum of operations
//! needed to transfrom one string into the other. An Operation is defined as an insertion, deletion,
//! substitution of a single character or transposition of two adjacent characters.
//!
//! Opposed to the [`Optimal String Alignment`] it allows editing a substring more than once.
//! An example where this leads to different results are the strings `CA` and `ÀBC`
//!
//! ```
//! use rapidfuzz::distance::damerau_levenshtein;
//! use rapidfuzz::distance::osa;
//!
//! assert_eq!(2, damerau_levenshtein::distance("CA".chars(), "ABC".chars(), None, None));
//! assert_eq!(3, osa::distance("CA".chars(), "ABC".chars(), None, None));
//! ```
//!
//! It does respect triangle inequality, and is thus a metric distance.
//!
//! [`Levenshtein`]: ../levenshtein/index.html
//! [`Optimal String Alignment`]: ../osa/index.html
//!
//! # Performance
//!
//! The implementation has a runtime complexity of `O(N*M)` and a memory usage of `O(N)`.
//! It's based on the paper
//! `Linear space string correction algorithm using the Damerau-Levenshtein distance`
//! from Chunchun Zhao and Sartaj Sahni
//!
//! ![benchmark results](https://raw.githubusercontent.com/maxbachmann/rapidfuzz-rs/main/doc/bench/damerau_levenshtein.svg)
//!

use crate::details::common::{norm_sim_to_norm_dist, remove_common_affix, HashableChar};
use crate::details::distance::{
    build_distance_metric_funcs, build_normalized_metric_funcs, CachedDistanceMetricUsize,
    CachedNormalizedDistanceMetricUsize,
};
use crate::details::growing_hashmap::HybridGrowingHashmap;
use std::cmp::{max, min};
use std::mem;

#[derive(Clone, Copy, PartialEq, Eq)]
struct RowId {
    val: isize,
}

impl Default for RowId {
    fn default() -> Self {
        RowId { val: -1 }
    }
}

/// based on the paper
/// "Linear space string correction algorithm using the Damerau-Levenshtein distance"
/// from Chunchun Zhao and Sartaj Sahni
///
/// todo in c++ this is templated on an integer type which reduced
/// memory usage depending on string lengths
fn damerau_damerau_levenshtein_distance_zhao<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> usize
where
    Iter1: Iterator<Item = Elem1>,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
{
    let max_val = max(len1, len2) as isize + 1;

    let mut last_row_id = HybridGrowingHashmap::<RowId> {
        map_unsigned: Default::default(),
        map_signed: Default::default(),
        extended_ascii: [Default::default(); 256],
    };
    let size = len2 + 2;
    let mut fr = vec![max_val; size];
    let mut r1 = vec![max_val; size];
    let mut r: Vec<isize> = (max_val..max_val + 1)
        .chain(0..(size - 1) as isize)
        .collect();

    for (i, ch1) in s1.enumerate().map(|(i, ch1)| (i + 1, ch1)) {
        mem::swap(&mut r, &mut r1);
        let mut last_col_id: isize = -1;
        let mut last_i2l1 = r[1];
        r[1] = i as isize;
        let mut t = max_val;

        for (j, ch2) in s2.clone().enumerate().map(|(j, ch2)| (j + 1, ch2)) {
            let diag = r1[j] + (ch1 != ch2) as isize;
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

    let dist = r[len2 + 1] as usize;

    if dist < score_cutoff {
        dist
    } else {
        score_cutoff + 1
    }
}

fn damerau_damerau_levenshtein_distance_impl<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> usize
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
{
    if score_cutoff < len1.abs_diff(len2) {
        return score_cutoff + 1;
    }

    let affix = remove_common_affix(s1, len1, s2, len2);
    damerau_damerau_levenshtein_distance_zhao(
        affix.s1,
        affix.len1,
        affix.s2,
        affix.len2,
        score_cutoff,
    )
}

pub(crate) struct DamerauLevenshtein {}

impl DamerauLevenshtein {
    build_distance_metric_funcs!(DamerauLevenshtein, usize, 0, usize::MAX);

    fn maximum(len1: usize, len2: usize) -> usize {
        max(len1, len2)
    }

    fn distance<Iter1, Iter2, Elem1, Elem2>(
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        damerau_damerau_levenshtein_distance_impl(s1, len1, s2, len2, score_cutoff)
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
/// assert_eq!(2, damerau_levenshtein::distance("CA".chars(), "ABC".chars(), None, None));
/// ```
pub fn distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> usize
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
    DamerauLevenshtein::distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(usize::MAX),
        score_hint.into().unwrap_or(usize::MAX),
    )
}

/// Damerau-Levenshtein similarity in the range [max, 0.0]
///
/// This is calculated as `max(len1, len2) - `[`distance`].
///
pub fn similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> usize
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
    DamerauLevenshtein::similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(0),
        score_hint.into().unwrap_or(0),
    )
}

/// Normalized Damerau-Levenshtein distance in the range [1.0, 0.0]
///
/// This is calculated as [`distance`]` / max(len1, len2)`.
///
pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> f64
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
    DamerauLevenshtein::normalized_distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(1.0),
        score_hint.into().unwrap_or(1.0),
    )
}

/// Normalized Damerau-Levenshtein similarity in the range [0.0, 1.0]
///
/// This is calculated as `1.0 - `[`normalized_distance`].
///
pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> f64
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
    DamerauLevenshtein::normalized_similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(0.0),
        score_hint.into().unwrap_or(0.0),
    )
}

/// Cached damerau levenshtein distance for `One x Many` comparisions
///
/// # Examples
///
/// ```
/// use rapidfuzz::distance::damerau_levenshtein;
///
/// let scorer = damerau_levenshtein::CachedDamerauLevenshtein::new("CA".chars());
/// assert_eq!(2, scorer.distance("ABC".chars(), None, None));
/// ```
pub struct CachedDamerauLevenshtein<Elem1> {
    s1: Vec<Elem1>,
}

impl<Elem1> CachedDistanceMetricUsize<Elem1> for CachedDamerauLevenshtein<Elem1>
where
    Elem1: HashableChar + Clone,
{
    fn maximum(&self, len2: usize) -> usize {
        max(self.s1.len(), len2)
    }

    fn _distance<Iter2, Elem2>(
        &self,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        damerau_damerau_levenshtein_distance_impl(
            self.s1.iter().copied(),
            self.s1.len(),
            s2,
            len2,
            score_cutoff,
        )
    }
}

impl<Elem1> CachedDamerauLevenshtein<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();
        CachedDamerauLevenshtein { s1 }
    }

    pub fn normalized_distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> f64
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        let len2 = s2_iter.clone().count();
        self._normalized_distance(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(1.0),
            score_hint.into().unwrap_or(1.0),
        )
    }

    pub fn normalized_similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> f64
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        let len2 = s2_iter.clone().count();
        self._normalized_similarity(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(0.0),
            score_hint.into().unwrap_or(0.0),
        )
    }

    /// test
    pub fn distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> usize
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        let s2_iter = s2.into_iter();
        let len2 = s2_iter.clone().count();
        self._distance(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(usize::MAX),
            score_hint.into().unwrap_or(usize::MAX),
        )
    }

    pub fn similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> usize
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        let s2_iter = s2.into_iter();
        let len2 = s2_iter.clone().count();
        self._similarity(
            s2_iter,
            len2,
            score_cutoff.into().unwrap_or(0),
            score_hint.into().unwrap_or(0),
        )
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
            if ($x - $y).abs() > $d {
                panic!();
            }
        };
    }

    fn _test_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> usize
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = distance(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = distance(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );

        let scorer1 = CachedDamerauLevenshtein::new(s1.clone());
        let res3 = scorer1.distance(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = CachedDamerauLevenshtein::new(s2.clone());
        let res4 = scorer2.distance(s1.clone(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn _test_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> usize
    where
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let res1 = _test_distance(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = _test_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    fn _test_normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> f64
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = normalized_similarity(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = normalized_similarity(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let scorer1 = CachedDamerauLevenshtein::new(s1.clone());
        let res3 =
            scorer1.normalized_similarity(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = CachedDamerauLevenshtein::new(s2.clone());
        let res4 = scorer2.normalized_similarity(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn _test_normalized_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> f64
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = _test_normalized_similarity(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = _test_normalized_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    /// levenshtein calculates empty sequence
    #[test]
    fn damerau_levenshtein_empty() {
        assert_eq!(0, _test_distance_ascii(EMPTY, EMPTY, None, None));
        assert_eq!(4, _test_distance_ascii(TEST, EMPTY, None, None));
        assert_eq!(4, _test_distance_ascii(EMPTY, TEST, None, None));
    }

    /// levenshtein calculates correct distances
    #[test]
    fn damerau_levenshtein_simple() {
        assert_eq!(0, _test_distance_ascii(TEST, TEST, None, None));
        assert_eq!(1, _test_distance_ascii(TEST, NO_SUFFIX, None, None));
        assert_eq!(1, _test_distance_ascii(TEST, NO_SUFFIX2, None, None));
        assert_eq!(1, _test_distance_ascii(SWAPPED1, SWAPPED2, None, None));
        assert_eq!(4, _test_distance_ascii(TEST, REPLACE_ALL, None, None));
        assert_eq!(2, _test_distance_ascii("CA", "ABC", None, None));

        assert_delta!(
            1.0,
            _test_normalized_similarity_ascii(TEST, TEST, None, None),
            0.0001
        );
        assert_delta!(
            0.75,
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX, None, None),
            0.0001
        );
        assert_delta!(
            0.75,
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX2, None, None),
            0.0001
        );
        assert_delta!(
            0.75,
            _test_normalized_similarity_ascii(SWAPPED1, SWAPPED2, None, None),
            0.0001
        );
        assert_delta!(
            0.0,
            _test_normalized_similarity_ascii(TEST, REPLACE_ALL, None, None),
            0.0001
        );
    }
}
