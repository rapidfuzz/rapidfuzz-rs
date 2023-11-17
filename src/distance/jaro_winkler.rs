use crate::details::common::{norm_sim_to_norm_dist, HashableChar, UnrefIterator};
use crate::details::distance::{
    build_cached_normalized_metric_funcs, build_cached_similarity_metric_funcs,
    build_normalized_metric_funcs, build_similarity_metric_funcs,
    less_than_score_cutoff_similarity,
};
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

pub(crate) struct JaroWinkler {}

impl JaroWinkler {
    build_similarity_metric_funcs!(JaroWinkler, f64, 0.0, 1.0, prefix_weight: Option<f64>);

    fn maximum<Iter1, Iter2, Elem1, Elem2>(
        _s1: Iter1,
        _len1: usize,
        _s2: Iter2,
        _len2: usize,
        _prefix_weight: Option<f64>,
    ) -> f64
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
        Elem1: PartialEq<Elem2>,
        Elem2: PartialEq<Elem1>,
    {
        1.0
    }

    pub(crate) fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        prefix_weight: Option<f64>,
        score_cutoff: f64,
        _score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + Clone + DoubleEndedIterator,
        Iter2: Iterator<Item = Elem2> + Clone + DoubleEndedIterator,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        jaro_winkler_similarity_without_pm(
            s1,
            len1,
            s2,
            len2,
            score_cutoff,
            prefix_weight.unwrap_or(0.1),
        )
    }
}

pub fn distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: Option<f64>,
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
    JaroWinkler::distance(s1, s2, prefix_weight, score_cutoff, score_hint)
}

pub fn similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: Option<f64>,
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
    JaroWinkler::similarity(s1, s2, prefix_weight, score_cutoff, score_hint)
}

pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: Option<f64>,
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
    JaroWinkler::normalized_distance(s1, s2, prefix_weight, score_cutoff, score_hint)
}

pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    prefix_weight: Option<f64>,
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
    JaroWinkler::normalized_similarity(s1, s2, prefix_weight, score_cutoff, score_hint)
}

pub struct CachedJaroWinkler<Elem1>
where
    Elem1: HashableChar + Clone,
{
    s1: Vec<Elem1>,
    pm: BlockPatternMatchVector,
    prefix_weight: f64,
}

impl<Elem1> CachedJaroWinkler<Elem1>
where
    Elem1: HashableChar + Clone,
{
    build_cached_similarity_metric_funcs!(CachedJaro, f64, 0.0, 1.0);

    pub fn new<Iter1>(s1: Iter1, prefix_weight: Option<f64>) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();

        let mut pm = BlockPatternMatchVector::new(s1.len());
        pm.insert(s1_iter);

        CachedJaroWinkler {
            s1,
            pm,
            prefix_weight: prefix_weight.unwrap_or(0.1),
        }
    }

    fn maximum<Iter2, Elem2>(&self, _s2: Iter2, _len2: usize) -> f64
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
    {
        1.0
    }

    fn _similarity<Iter2, Elem2>(
        &self,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        _score_hint: f64,
    ) -> f64
    where
        Iter2: Iterator<Item = Elem2> + Clone + DoubleEndedIterator,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        jaro_winkler_similarity_with_pm(
            &self.pm,
            UnrefIterator {
                seq: self.s1.iter(),
            },
            self.s1.len(),
            s2,
            len2,
            self.prefix_weight,
            score_cutoff,
        )
    }
}
