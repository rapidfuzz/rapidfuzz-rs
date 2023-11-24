use crate::details::common::{find_common_prefix, HashableChar};
use crate::details::distance::{NormalizedMetricUsize2, SimilarityMetricUsize};

struct Prefix;

impl SimilarityMetricUsize for Prefix {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        _len1: usize,
        s2: Iter2,
        _len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + Clone,
        Iter2: Iterator<Item = Elem2> + Clone,
        Elem1: PartialEq<Elem2> + HashableChar,
        Elem2: PartialEq<Elem1> + HashableChar,
    {
        let dist = find_common_prefix(s1, s2);
        if dist >= score_cutoff {
            dist
        } else {
            0
        }
    }
}

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
    Prefix {}._distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(usize::MAX),
        score_hint.into().unwrap_or(usize::MAX),
    )
}

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
    Prefix {}._similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(0),
        score_hint.into().unwrap_or(0),
    )
}

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
    Prefix {}._normalized_distance(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(1.0),
        score_hint.into().unwrap_or(1.0),
    )
}

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
    Prefix {}._normalized_similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into().unwrap_or(0.0),
        score_hint.into().unwrap_or(0.0),
    )
}

pub struct CachedPrefix<Elem1> {
    s1: Vec<Elem1>,
}

impl<Elem1> CachedPrefix<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1.into_iter();
        Self {
            s1: s1_iter.collect(),
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
        normalized_distance(self.s1.iter().copied(), s2, score_cutoff, score_hint)
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
        normalized_similarity(self.s1.iter().copied(), s2, score_cutoff, score_hint)
    }

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
        distance(self.s1.iter().copied(), s2, score_cutoff, score_hint)
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
        similarity(self.s1.iter().copied(), s2, score_cutoff, score_hint)
    }
}

// todo _or versions
