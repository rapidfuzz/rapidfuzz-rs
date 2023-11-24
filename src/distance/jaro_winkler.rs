use crate::details::common::HashableChar;
use crate::details::distance::{NormalizedMetricf64, SimilarityMetricf64};
use crate::details::pattern_match_vector::BlockPatternMatchVector;

use crate::distance::jaro::{jaro_similarity_with_pm, jaro_similarity_without_pm};

fn jaro_winkler_similarity_without_pm<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    prefix_weight: f64,
    score_cutoff: f64,
) -> f64
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
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

    let mut sim = jaro_similarity_without_pm(s1, len1, s2, len2, jaro_score_cutoff);
    if sim > 0.7 {
        sim += prefix as f64 * prefix_weight / (1.0 - sim);
    }

    if sim >= score_cutoff {
        sim
    } else {
        0.0
    }
}

fn jaro_winkler_similarity_with_pm<Iter1, Iter2, Elem1, Elem2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    prefix_weight: f64,
    score_cutoff: f64,
) -> f64
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
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

    let mut sim = jaro_similarity_with_pm(pm, s1, len1, s2, len2, jaro_score_cutoff);
    if sim > 0.7 {
        sim += prefix as f64 * prefix_weight / (1.0 - sim);
    }

    if sim >= score_cutoff {
        sim
    } else {
        0.0
    }
}

pub(crate) struct JaroWinkler {
    prefix_weight: Option<f64>,
}

impl SimilarityMetricf64 for JaroWinkler {
    fn maximum(&self, _len1: usize, _len2: usize) -> f64 {
        1.0
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        _score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        jaro_winkler_similarity_without_pm(
            s1,
            len1,
            s2,
            len2,
            score_cutoff,
            self.prefix_weight.unwrap_or(0.1),
        )
    }
}

pub fn distance<Iter1, Iter2, Elem1, Elem2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
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
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    JaroWinkler {
        prefix_weight: prefix_weight.into(),
    }
    ._distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(1.0),
        score_hint.into().unwrap_or(1.0),
    )
}

pub fn similarity<Iter1, Iter2, Elem1, Elem2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
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
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    JaroWinkler {
        prefix_weight: prefix_weight.into(),
    }
    ._similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(0.0),
        score_hint.into().unwrap_or(0.0),
    )
}

pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
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
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    JaroWinkler {
        prefix_weight: prefix_weight.into(),
    }
    ._normalized_distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(1.0),
        score_hint.into().unwrap_or(1.0),
    )
}

pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2, PrefixWeight, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: PrefixWeight,
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
    PrefixWeight: Into<Option<f64>>,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    JaroWinkler {
        prefix_weight: prefix_weight.into(),
    }
    ._normalized_similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(0.0),
        score_hint.into().unwrap_or(0.0),
    )
}

pub struct CachedJaroWinkler<Elem1> {
    s1: Vec<Elem1>,
    pm: BlockPatternMatchVector,
    prefix_weight: f64,
}

impl<CharT> SimilarityMetricf64 for CachedJaroWinkler<CharT> {
    fn maximum(&self, _len1: usize, _len2: usize) -> f64 {
        1.0
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        _score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        jaro_winkler_similarity_with_pm(
            &self.pm,
            s1,
            len1,
            s2,
            len2,
            self.prefix_weight,
            score_cutoff,
        )
    }
}

impl<Elem1> CachedJaroWinkler<Elem1>
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
        self._normalized_distance(
            self.s1.iter().copied(),
            self.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
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
        self._normalized_similarity(
            self.s1.iter().copied(),
            self.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into().unwrap_or(0.0),
            score_hint.into().unwrap_or(0.0),
        )
    }

    pub fn distance<Iter2, Elem2, ScoreCutoff, ScoreHint>(
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
        self._distance(
            self.s1.iter().copied(),
            self.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into().unwrap_or(1.0),
            score_hint.into().unwrap_or(1.0),
        )
    }

    pub fn similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
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
        self._similarity(
            self.s1.iter().copied(),
            self.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff.into().unwrap_or(0.0),
            score_hint.into().unwrap_or(0.0),
        )
    }
}
