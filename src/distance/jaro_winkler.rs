//! Jaro-Winkler similarity
//!
//! The Jaro-Winkler similarity extends the [`Jaro`] similarity to provide additional
//! sensitivity to matching prefixes. It introduces a scaling mechanism that boosts
//! the similarity score for strings with common prefixes.
//!
//! [`Jaro`]: ../jaro/index.html
//!
//! # Performance
//!
//! The implementation has a runtime complexity of `O([N/64]*M)` and a memory usage of `O(N)`.
//!
//! ![benchmark results](https://raw.githubusercontent.com/maxbachmann/rapidfuzz-rs/main/rapidfuzz-benches/results/jaro_winkler.svg)
//!

use crate::details::distance::Metricf64;
use crate::details::pattern_match_vector::BlockPatternMatchVector;
use crate::HashableChar;

use crate::distance::jaro;

fn similarity_without_pm<Iter1, Iter2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    prefix_weight: f64,
    score_cutoff: f64,
) -> Option<f64>
where
    Iter1: Iterator + Clone,
    Iter2: Iterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    let prefix = s1
        .clone()
        .zip(s2.clone())
        .take(4)
        .take_while(|(ch1, ch2)| ch1 == ch2)
        .count();

    let mut jaro_score_cutoff = score_cutoff;
    if jaro_score_cutoff > 0.7 {
        let prefix_sim = prefix as f64 * prefix_weight;
        jaro_score_cutoff = if prefix_sim >= 1.0 {
            0.7
        } else {
            0.7_f64.max((prefix_sim - jaro_score_cutoff) / (prefix_sim - 1.0))
        }
    }

    let mut sim = jaro::similarity_without_pm(s1, len1, s2, len2, jaro_score_cutoff)?;
    if sim > 0.7 {
        sim += prefix as f64 * prefix_weight * (1.0 - sim);
    }

    if sim >= score_cutoff {
        Some(sim)
    } else {
        None
    }
}

fn similarity_with_pm<Iter1, Iter2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    prefix_weight: f64,
    score_cutoff: f64,
) -> Option<f64>
where
    Iter1: Iterator + Clone,
    Iter2: Iterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    let prefix = s1
        .clone()
        .zip(s2.clone())
        .take(4)
        .take_while(|(ch1, ch2)| ch1 == ch2)
        .count();

    let mut jaro_score_cutoff = score_cutoff;
    if jaro_score_cutoff > 0.7 {
        let prefix_sim = prefix as f64 * prefix_weight;
        jaro_score_cutoff = if prefix_sim >= 1.0 {
            0.7
        } else {
            0.7_f64.max((prefix_sim - jaro_score_cutoff) / (prefix_sim - 1.0))
        }
    }

    let mut sim = jaro::similarity_with_pm(pm, s1, len1, s2, len2, jaro_score_cutoff)?;
    if sim > 0.7 {
        sim += prefix as f64 * prefix_weight * (1.0 - sim);
    }

    if sim >= score_cutoff {
        Some(sim)
    } else {
        None
    }
}

pub(crate) struct IndividualComparator {
    prefix_weight: Option<f64>,
}

impl Metricf64 for IndividualComparator {
    fn maximum(&self, _len1: usize, _len2: usize) -> f64 {
        1.0
    }

    fn _similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        _score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        similarity_without_pm(
            s1,
            len1,
            s2,
            len2,
            self.prefix_weight.unwrap_or(0.1),
            score_cutoff.unwrap_or(0.0),
        )
    }
}

/// Jaro-Winkler distance in the range [0.0, 1.0].
///
/// This is calculated as `1.0 - `[`similarity`].
///
pub fn distance<Iter1, Iter2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    IndividualComparator {
        prefix_weight: prefix_weight.into(),
    }
    ._distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

/// Jaro-Winkler similarity in the range [1.0, 0.0].
pub fn similarity<Iter1, Iter2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    IndividualComparator {
        prefix_weight: prefix_weight.into(),
    }
    ._similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

/// Normalized Jaro-Winkler distance in the range [0.0, 1.0].
///
/// This behaves the same as `distance`, since the Jaro-Winkler similarity is always
/// normalized
///
pub fn normalized_distance<Iter1, Iter2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    IndividualComparator {
        prefix_weight: prefix_weight.into(),
    }
    ._normalized_distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

/// Normalized Jaro-Winkler similarity in the range [1.0, 0.0].
///
/// This behaves the same as `similarity`, since the Jaro-Winkler similarity is always
/// normalized
///
pub fn normalized_similarity<Iter1, Iter2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    IndividualComparator {
        prefix_weight: prefix_weight.into(),
    }
    ._normalized_similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

/// `One x Many` comparisons using the Jaro-Winkler similarity
pub struct BatchComparator<Elem1> {
    s1: Vec<Elem1>,
    pm: BlockPatternMatchVector,
    prefix_weight: f64,
}

impl<CharT> Metricf64 for BatchComparator<CharT> {
    fn maximum(&self, _len1: usize, _len2: usize) -> f64 {
        1.0
    }

    fn _similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        _score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        similarity_with_pm(
            &self.pm,
            s1,
            len1,
            s2,
            len2,
            self.prefix_weight,
            score_cutoff.unwrap_or(0.0),
        )
    }
}

impl<Elem1> BatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1_: Iter1, prefix_weight: Option<f64>) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1_.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();

        let mut pm = BlockPatternMatchVector::new(s1.len());
        pm.insert(s1_iter);

        Self {
            s1,
            pm,
            prefix_weight: prefix_weight.unwrap_or(0.1),
        }
    }

    /// Normalized distance calculated similar to [`normalized_distance`]
    pub fn normalized_distance<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        self._normalized_distance(
            self.s1.iter().copied(),
            self.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into(),
            score_hint.into(),
        )
    }

    /// Normalized similarity calculated similar to [`normalized_similarity`]
    pub fn normalized_similarity<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        self._normalized_similarity(
            self.s1.iter().copied(),
            self.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into(),
            score_hint.into(),
        )
    }

    /// Distance calculated similar to [`distance`]
    pub fn distance<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        self._distance(
            self.s1.iter().copied(),
            self.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into(),
            score_hint.into(),
        )
    }

    /// Similarity calculated similar to [`similarity`]
    pub fn similarity<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
    {
        let s2_iter = s2.into_iter();
        self._similarity(
            self.s1.iter().copied(),
            self.s1.len(),
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

    fn _test_distance<Iter1, Iter2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = distance(
            s1.clone(),
            s2.clone(),
            None,
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = distance(
            s2.clone(),
            s1.clone(),
            None,
            score_cutoff.clone(),
            score_hint.clone(),
        );

        let scorer1 = BatchComparator::new(s1.clone(), None);
        let res3 = scorer1.distance(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone(), None);
        let res4 = scorer2.distance(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn _test_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = _test_distance(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = _test_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    fn _test_similarity<Iter1, Iter2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = similarity(
            s1.clone(),
            s2.clone(),
            None,
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = similarity(
            s2.clone(),
            s1.clone(),
            None,
            score_cutoff.clone(),
            score_hint.clone(),
        );

        let scorer1 = BatchComparator::new(s1.clone(), None);
        let res3 = scorer1.similarity(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone(), None);
        let res4 = scorer2.similarity(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn _test_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = _test_similarity(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = _test_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    #[test]
    fn test_no_cutoff() {
        assert_delta!(
            Some(0.455556),
            _test_similarity_ascii("james", "robert", None, None),
            0.0001
        );
        assert_delta!(
            Some(1.0 - 0.455556),
            _test_distance_ascii("james", "robert", None, None),
            0.0001
        );
    }

    #[test]
    fn test_flag_chars() {
        let names = [
            "james",
            "robert",
            "john",
            "michael",
            "william",
            "david",
            "joseph",
            "thomas",
            "charles",
            "mary",
            "patricia",
            "jennifer",
            "linda",
            "elizabeth",
            "barbara",
            "susan",
            "jessica",
            "sarah",
            "karen",
            "",
            "aaaaaaaa",
            "aabaaab",
        ];

        let score_cutoffs = [0.0]; //, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1];

        let scores = [
            1.0, 0.455556, 0.483333, 0.561905, 0.0, 0.466667, 0.588889, 0.577778, 0.67619,
            0.483333, 0.441667, 0.55, 0.0, 0.374074, 0.447619, 0.0, 0.67619, 0.466667, 0.6, 0.0,
            0.441667, 0.447619, 0.455556, 1.0, 0.472222, 0.436508, 0.0, 0.0, 0.555556, 0.444444,
            0.373016, 0.472222, 0.361111, 0.527778, 0.0, 0.5, 0.531746, 0.0, 0.436508, 0.455556,
            0.577778, 0.0, 0.0, 0.436508, 0.483333, 0.472222, 1.0, 0.464286, 0.0, 0.0, 0.611111,
            0.444444, 0.464286, 0.0, 0.0, 0.583333, 0.483333, 0.0, 0.0, 0.483333, 0.464286, 0.0,
            0.483333, 0.0, 0.0, 0.0, 0.561905, 0.436508, 0.464286, 1.0, 0.52381, 0.447619,
            0.373016, 0.539683, 0.742857, 0.464286, 0.490079, 0.511905, 0.561905, 0.587302,
            0.428571, 0.447619, 0.428571, 0.395238, 0.447619, 0.0, 0.422619, 0.428571, 0.0, 0.0,
            0.0, 0.52381, 1.0, 0.447619, 0.0, 0.436508, 0.428571, 0.0, 0.60119, 0.422619, 0.565079,
            0.47619, 0.428571, 0.447619, 0.52381, 0.447619, 0.0, 0.0, 0.422619, 0.428571, 0.466667,
            0.0, 0.0, 0.447619, 0.447619, 1.0, 0.0, 0.0, 0.447619, 0.483333, 0.55, 0.441667,
            0.466667, 0.374074, 0.447619, 0.0, 0.447619, 0.466667, 0.466667, 0.0, 0.441667,
            0.447619, 0.588889, 0.555556, 0.611111, 0.373016, 0.0, 0.0, 1.0, 0.444444, 0.436508,
            0.0, 0.0, 0.527778, 0.0, 0.518519, 0.0, 0.455556, 0.531746, 0.577778, 0.455556, 0.0,
            0.0, 0.0, 0.577778, 0.444444, 0.444444, 0.539683, 0.436508, 0.0, 0.444444, 1.0,
            0.642857, 0.0, 0.361111, 0.0, 0.455556, 0.425926, 0.436508, 0.455556, 0.373016,
            0.455556, 0.0, 0.0, 0.430556, 0.436508, 0.67619, 0.373016, 0.464286, 0.742857,
            0.428571, 0.447619, 0.436508, 0.642857, 1.0, 0.595238, 0.511905, 0.422619, 0.447619,
            0.47619, 0.52381, 0.447619, 0.0, 0.561905, 0.67619, 0.0, 0.422619, 0.428571, 0.483333,
            0.472222, 0.0, 0.464286, 0.0, 0.483333, 0.0, 0.0, 0.595238, 1.0, 0.583333, 0.0, 0.0,
            0.453704, 0.595238, 0.0, 0.0, 0.633333, 0.633333, 0.0, 0.458333, 0.464286, 0.441667,
            0.361111, 0.0, 0.490079, 0.60119, 0.55, 0.0, 0.361111, 0.511905, 0.583333, 1.0,
            0.416667, 0.383333, 0.324074, 0.60119, 0.441667, 0.60119, 0.55, 0.55, 0.0, 0.5,
            0.511905, 0.55, 0.527778, 0.583333, 0.511905, 0.422619, 0.441667, 0.527778, 0.0,
            0.422619, 0.0, 0.416667, 1.0, 0.383333, 0.569444, 0.422619, 0.441667, 0.60119, 0.0,
            0.55, 0.0, 0.0, 0.0, 0.0, 0.0, 0.483333, 0.561905, 0.565079, 0.466667, 0.0, 0.455556,
            0.447619, 0.0, 0.383333, 0.383333, 1.0, 0.644444, 0.447619, 0.466667, 0.447619,
            0.466667, 0.0, 0.0, 0.441667, 0.447619, 0.374074, 0.5, 0.0, 0.587302, 0.47619,
            0.374074, 0.518519, 0.425926, 0.47619, 0.453704, 0.324074, 0.569444, 0.644444, 1.0,
            0.502646, 0.437037, 0.587302, 0.437037, 0.374074, 0.0, 0.412037, 0.502646, 0.447619,
            0.531746, 0.0, 0.428571, 0.428571, 0.447619, 0.0, 0.436508, 0.52381, 0.595238, 0.60119,
            0.422619, 0.447619, 0.502646, 1.0, 0.447619, 0.428571, 0.67619, 0.561905, 0.0, 0.60119,
            0.630952, 0.0, 0.0, 0.483333, 0.447619, 0.447619, 0.0, 0.455556, 0.455556, 0.447619,
            0.0, 0.441667, 0.441667, 0.466667, 0.437037, 0.447619, 1.0, 0.561905, 0.6, 0.466667,
            0.0, 0.441667, 0.447619, 0.67619, 0.436508, 0.464286, 0.428571, 0.52381, 0.447619,
            0.531746, 0.373016, 0.0, 0.0, 0.60119, 0.60119, 0.447619, 0.587302, 0.428571, 0.561905,
            1.0, 0.447619, 0.447619, 0.0, 0.422619, 0.428571, 0.466667, 0.455556, 0.0, 0.395238,
            0.447619, 0.466667, 0.577778, 0.455556, 0.561905, 0.633333, 0.55, 0.0, 0.466667,
            0.437037, 0.67619, 0.6, 0.447619, 1.0, 0.6, 0.0, 0.55, 0.561905, 0.6, 0.577778,
            0.483333, 0.447619, 0.0, 0.466667, 0.455556, 0.0, 0.67619, 0.633333, 0.55, 0.55, 0.0,
            0.374074, 0.561905, 0.466667, 0.447619, 0.6, 1.0, 0.0, 0.441667, 0.447619, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 0.441667, 0.0, 0.0, 0.422619, 0.422619, 0.441667, 0.0, 0.430556,
            0.422619, 0.458333, 0.5, 0.0, 0.441667, 0.412037, 0.60119, 0.441667, 0.422619, 0.55,
            0.441667, 0.0, 1.0, 0.82381, 0.447619, 0.436508, 0.0, 0.428571, 0.428571, 0.447619,
            0.0, 0.436508, 0.428571, 0.464286, 0.511905, 0.0, 0.447619, 0.502646, 0.630952,
            0.447619, 0.428571, 0.561905, 0.447619, 0.0, 0.82381, 1.0,
        ];
        for score_cutoff in score_cutoffs {
            for (i, name1) in names.iter().enumerate() {
                for (j, name2) in names.iter().enumerate() {
                    let score = scores[i * names.len() + j];
                    let expected_sim = if score_cutoff <= score {
                        Some(score)
                    } else {
                        None
                    };
                    let expected_dist = expected_sim.map(|s| 1.0 - s);

                    let sim = _test_similarity_ascii(name1, name2, score_cutoff, None);
                    let dist = _test_distance_ascii(name1, name2, 1.0 - score_cutoff, None);
                    assert_delta!(expected_sim, sim, 0.0001);
                    assert_delta!(expected_dist, dist, 0.0001);
                }
            }
        }
    }

    #[test]
    fn unicode() {
        assert_delta!(
            Some(0.375),
            _test_distance("Иванко".chars(), "Петрунко".chars(), None, None),
            0.0001
        );
    }
}
