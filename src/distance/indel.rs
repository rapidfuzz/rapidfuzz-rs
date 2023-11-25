use crate::details::common::HashableChar;
use crate::details::distance::MetricUsize;
use crate::details::pattern_match_vector::BlockPatternMatchVector;
use crate::distance::lcs_seq::{lcs_seq_similarity_with_pm, CachedLcsSeq, LcsSeq};

pub(crate) struct Indel;

impl MetricUsize for Indel {
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
            LcsSeq {}._similarity(s1, len1, s2, len2, Some(lcs_cutoff), Some(lcs_hint))?;
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
    Indel {}._distance(
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
    Indel {}._similarity(
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
    Indel {}._normalized_distance(
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
    Indel {}._normalized_similarity(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        score_cutoff.into(),
        score_hint.into(),
    )
}

pub(crate) fn indel_distance_with_pm<Iter1, Iter2, Elem1, Elem2>(
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

    let lcs_sim = lcs_seq_similarity_with_pm(pm, s1, len1, s2, len2, lcs_cutoff)?;
    let dist = maximum - 2 * lcs_sim;
    if dist <= score_cutoff {
        Some(dist)
    } else {
        None
    }
}

pub struct CachedIndel<Elem1>
where
    Elem1: HashableChar + Clone,
{
    scorer: CachedLcsSeq<Elem1>,
}

impl<CharT> MetricUsize for CachedIndel<CharT>
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

impl<Elem1> CachedIndel<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        Self {
            scorer: CachedLcsSeq::new(s1),
        }
    }

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
