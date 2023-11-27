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
//! ![benchmark results](https://raw.githubusercontent.com/maxbachmann/rapidfuzz-rs/main/rapidfuzz-benches/results/indel.svg)
//!

use crate::details::distance::MetricUsize;
use crate::details::pattern_match_vector::BlockPatternMatchVector;
use crate::distance::lcs_seq;
use crate::HashableChar;

pub(crate) struct IndividualComparator;

impl MetricUsize for IndividualComparator {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1 + len2
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        score_hint_: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
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
        )?;
        let dist = maximum - 2 * lcs_sim;
        if dist <= score_cutoff {
            Some(dist)
        } else {
            None
        }
    }
}

pub fn distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<usize>
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
    IndividualComparator {}._distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

pub fn similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<usize>
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
    IndividualComparator {}._similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
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
    IndividualComparator {}._normalized_distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
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
    IndividualComparator {}._normalized_similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

pub(crate) fn distance_with_pm<Iter1, Iter2, Elem1, Elem2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    let maximum = len1 + len2;
    let lcs_cutoff = if maximum / 2 >= score_cutoff {
        maximum / 2 - score_cutoff
    } else {
        0
    };

    let lcs_sim = lcs_seq::similarity_with_pm(pm, s1, len1, s2, len2, lcs_cutoff)?;
    let dist = maximum - 2 * lcs_sim;
    if dist <= score_cutoff {
        Some(dist)
    } else {
        None
    }
}

pub struct BatchComparator<Elem1> {
    scorer: lcs_seq::BatchComparator<Elem1>,
}

impl<CharT> MetricUsize for BatchComparator<CharT>
where
    CharT: HashableChar + Clone,
{
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1 + len2
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        score_hint_: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
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
        let lcs_sim =
            self.scorer
                ._similarity(s1, len1, s2, len2, Some(lcs_cutoff), Some(lcs_hint))?;
        let dist = maximum - 2 * lcs_sim;
        if dist <= score_cutoff {
            Some(dist)
        } else {
            None
        }
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
    pub fn normalized_distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        self._normalized_distance(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into(),
            score_hint.into(),
        )
    }

    /// Normalized similarity calculated similar to [`normalized_similarity`]
    pub fn normalized_similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        self._normalized_similarity(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into(),
            score_hint.into(),
        )
    }

    /// Distance calculated similar to [`distance`]
    pub fn distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        let s2_iter = s2.into_iter();
        self._distance(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into(),
            score_hint.into(),
        )
    }

    /// Similarity calculated similar to [`similarity`]
    pub fn similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<usize>>,
        ScoreHint: Into<Option<usize>>,
    {
        let s2_iter = s2.into_iter();
        self._similarity(
            self.scorer.s1.iter().copied(),
            self.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into(),
            score_hint.into(),
        )
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

    fn test_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
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

        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.distance(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.distance(s1.clone(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn test_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
    where
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let res1 = test_distance(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    fn test_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
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
        let res1 = similarity(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = similarity(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );

        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.similarity(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.similarity(s1.clone(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn test_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
    where
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let res1 = test_similarity(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    fn test_normalized_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
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
        let res1 = normalized_distance(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = normalized_distance(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let scorer1 = BatchComparator::new(s1.clone());
        let res3 =
            scorer1.normalized_distance(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.normalized_distance(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn test_normalized_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = test_normalized_distance(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_normalized_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    fn test_normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
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
        let scorer1 = BatchComparator::new(s1.clone());
        let res3 =
            scorer1.normalized_similarity(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.normalized_similarity(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn test_normalized_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = test_normalized_similarity(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_normalized_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    #[test]
    fn similar() {
        assert_eq!(Some(0), test_distance_ascii("aaaa", "aaaa", None, None));
        assert_eq!(Some(8), test_similarity_ascii("aaaa", "aaaa", None, None));
        assert_delta!(
            Some(0.0),
            test_normalized_distance_ascii("aaaa", "aaaa", None, None),
            0.0001
        );
        assert_delta!(
            Some(1.0),
            test_normalized_similarity_ascii("aaaa", "aaaa", None, None),
            0.0001
        );
    }

    #[test]
    fn completely_different() {
        assert_eq!(Some(8), test_distance_ascii("aaaa", "bbbb", None, None));
        assert_eq!(Some(0), test_similarity_ascii("aaaa", "bbbb", None, None));
        assert_delta!(
            Some(1.0),
            test_normalized_distance_ascii("aaaa", "bbbb", None, None),
            0.0001
        );
        assert_delta!(
            Some(0.0),
            test_normalized_similarity_ascii("aaaa", "bbbb", None, None),
            0.0001
        );
    }

    #[test]
    fn test_mbleven() {
        let mut a = "South Korea";
        let mut b = "North Korea";

        assert_eq!(Some(4), test_distance_ascii(a, b, None, None));
        assert_eq!(Some(4), test_distance_ascii(a, b, 5, None));
        assert_eq!(Some(4), test_distance_ascii(a, b, 4, None));
        assert_eq!(None, test_distance_ascii(a, b, 3, None));
        assert_eq!(None, test_distance_ascii(a, b, 2, None));
        assert_eq!(None, test_distance_ascii(a, b, 1, None));
        assert_eq!(None, test_distance_ascii(a, b, 0, None));

        a = "aabc";
        b = "cccd";
        assert_eq!(Some(6), test_distance_ascii(a, b, None, None));
        assert_eq!(Some(6), test_distance_ascii(a, b, 6, None));
        assert_eq!(None, test_distance_ascii(a, b, 5, None));
        assert_eq!(None, test_distance_ascii(a, b, 4, None));
        assert_eq!(None, test_distance_ascii(a, b, 3, None));
        assert_eq!(None, test_distance_ascii(a, b, 2, None));
        assert_eq!(None, test_distance_ascii(a, b, 1, None));
        assert_eq!(None, test_distance_ascii(a, b, 0, None));
    }

    // this was an issue in the cached lcs implementation of rapidfuzz-cpp
    #[test]
    fn test_issue_unknown() {
        let s1 = "001";
        let s2 = "220";
        assert_delta!(
            Some(0.3333333),
            test_normalized_similarity_ascii(s1, s2, None, None),
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
        assert_eq!(Some(508), test_distance_ascii(s1, s2, None, None));
        assert_eq!(Some(508), test_distance_ascii(s1, s2, 508, None));
        assert_eq!(None, test_distance_ascii(s1, s2, 507, None));
        assert_eq!(Some(508), test_distance_ascii(s1, s2, usize::MAX, None));

        s1 = "bbbdbbmbbbbbbbbbBbfbbbbbbbbbbbbbbbbbbbrbbbbbrbbbbbdbnbbbjbhbbbbbbbbbhbbb\
            bbCbobbbxbbbbbkbbbAbxbbwbbbtbcbbbbebbiblbbbbqbbbbbbpbbbbbbubbbkbbDbbbhbkbC\
            bbgbbrbbbbbbbbbbbkbyvbbsbAbbbbz";
        s2 = "jaaagaaqyaaaanrCfwaaxaeahtaaaCzaaaspaaBkvaaaaqDaacndaaeolwiaaauaaaaaaamA";

        assert_eq!(Some(231), test_distance_ascii(s1, s2, None, None));
    }
}
