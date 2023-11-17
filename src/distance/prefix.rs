use crate::details::common::{
    find_common_prefix, norm_sim_to_norm_dist, HashableChar, UnrefIterator,
};
use crate::details::distance::{
    build_cached_normalized_metric_funcs, build_cached_similarity_metric_funcs,
    build_normalized_metric_funcs, build_similarity_metric_funcs,
    less_than_score_cutoff_similarity,
};

struct Prefix {}

impl Prefix {
    build_similarity_metric_funcs!(Prefix, usize, 0, usize::MAX);

    fn maximum<Iter1, Iter2, Elem1, Elem2>(
        _s1: Iter1,
        len1: usize,
        _s2: Iter2,
        len2: usize,
    ) -> usize
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
        Elem1: PartialEq<Elem2> + HashableChar,
        Elem2: PartialEq<Elem1> + HashableChar,
    {
        len1.max(len2)
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
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
    Prefix::distance(s1, s2, score_cutoff, score_hint)
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
    Prefix::similarity(s1, s2, score_cutoff, score_hint)
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
    Prefix::normalized_distance(s1, s2, score_cutoff, score_hint)
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
    Prefix::normalized_similarity(s1, s2, score_cutoff, score_hint)
}

pub struct CachedPrefix<Elem1>
where
    Elem1: HashableChar + Clone,
{
    s1: Vec<Elem1>,
}

impl<Elem1> CachedPrefix<Elem1>
where
    Elem1: HashableChar + Clone,
{
    build_cached_similarity_metric_funcs!(CachedPrefix, usize, 0, usize::MAX);

    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();
        CachedPrefix { s1 }
    }

    fn maximum<Iter2, Elem2>(&self, _s2: Iter2, len2: usize) -> usize
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
    {
        self.s1.len().max(len2)
    }

    fn _similarity<Iter2, Elem2>(
        &self,
        s2: Iter2,
        _len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter2: Iterator<Item = Elem2> + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let dist = find_common_prefix(
            UnrefIterator {
                seq: self.s1.iter(),
            },
            s2,
        );
        if dist >= score_cutoff {
            dist
        } else {
            0
        }
    }
}
