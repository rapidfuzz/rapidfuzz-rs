use crate::details::common::{norm_sim_to_norm_dist, HashableChar};
use crate::details::distance::{
    build_cached_distance_metric_funcs, build_cached_normalized_metric_funcs,
    build_distance_metric_funcs, build_normalized_metric_funcs,
};
use crate::details::pattern_match_vector::BlockPatternMatchVector;
use crate::distance::lcs_seq::{lcs_seq_similarity_with_pm, CachedLcsSeq, LcsSeq};

pub(crate) struct Indel {}

impl Indel {
    build_distance_metric_funcs!(Indel, usize, 0, usize::MAX);

    fn maximum(len1: usize, len2: usize) -> usize {
        len1 + len2
    }

    pub(crate) fn _distance<Iter1, Iter2, Elem1, Elem2>(
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + Clone + DoubleEndedIterator,
        Iter2: Iterator<Item = Elem2> + Clone + DoubleEndedIterator,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = Indel::maximum(len1, len2);
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
        let lcs_sim = LcsSeq::_similarity(s1, len1, s2, len2, lcs_cutoff, lcs_hint);
        let dist = maximum - 2 * lcs_sim;
        if dist <= score_cutoff {
            dist
        } else {
            score_cutoff + 1
        }
    }
}

pub fn distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: Option<usize>,
    score_hint: Option<usize>,
) -> usize
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    Indel::distance(s1, s2, score_cutoff, score_hint)
}

pub fn similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: Option<usize>,
    score_hint: Option<usize>,
) -> usize
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    Indel::similarity(s1, s2, score_cutoff, score_hint)
}

pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: Option<f64>,
    score_hint: Option<f64>,
) -> f64
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    Indel::normalized_distance(s1, s2, score_cutoff, score_hint)
}

pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: Option<f64>,
    score_hint: Option<f64>,
) -> f64
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
    <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    Indel::normalized_similarity(s1, s2, score_cutoff, score_hint)
}

pub(crate) fn indel_distance_with_pm<Iter1, Iter2, Elem1, Elem2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> usize
where
    Iter1: Iterator<Item = Elem1> + Clone + DoubleEndedIterator,
    Iter2: Iterator<Item = Elem2> + Clone + DoubleEndedIterator,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    let maximum = len1 + len2;
    let lcs_cutoff = if maximum / 2 >= score_cutoff {
        maximum / 2 - score_cutoff
    } else {
        0
    };

    let lcs_sim = lcs_seq_similarity_with_pm(pm, s1, len1, s2, len2, lcs_cutoff);
    let dist = maximum - 2 * lcs_sim;
    if dist <= score_cutoff {
        dist
    } else {
        score_cutoff + 1
    }
}

pub struct CachedIndel<Elem1>
where
    Elem1: HashableChar + Clone,
{
    scorer: CachedLcsSeq<Elem1>,
}

impl<Elem1> CachedIndel<Elem1>
where
    Elem1: HashableChar + Clone,
{
    build_cached_distance_metric_funcs!(CachedIndel, usize, 0, usize::MAX);

    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        CachedIndel {
            scorer: CachedLcsSeq::new(s1),
        }
    }

    fn maximum(&self, len2: usize) -> usize {
        self.scorer.s1.len() + len2
    }

    fn _distance<Iter2, Elem2>(
        &self,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        score_hint: usize,
    ) -> usize
    where
        Iter2: Iterator<Item = Elem2> + Clone + DoubleEndedIterator,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len2);
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
        let lcs_sim = self.scorer._similarity(s2, len2, lcs_cutoff, lcs_hint);
        let dist = maximum - 2 * lcs_sim;
        if dist <= score_cutoff {
            dist
        } else {
            score_cutoff + 1
        }
    }
}
