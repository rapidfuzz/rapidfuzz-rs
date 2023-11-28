//! Optimal String Alignment distance
//!
//! The Optimal String Alignment distance (OSA) measures the minimum number of operations required to
//! transform one string into another, considering four types of elementary edits:
//! `insertions`, `deletions`, `substitutions`, and `transpositions`
//!
//! # Differences from Damerau-Levenshtein distance
//!
//! While both the [`Damerau-Levenshtein`] and OSA distance include transpositions,
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
//! assert_eq!(Some(2), damerau_levenshtein::distance("CA".chars(), "ABC".chars(), None, None));
//! assert_eq!(3, osa::distance().compare("CA".chars(), "ABC".chars()));
//! ```
//!
//! The handling of transpositions in the OSA distance is simpler, which makes it computationally less intensive.
//!
//! ## Performance
//!
//! The implementation has a runtime complexity of `O([N/64]*M)` and a memory usage of `O(N)`.
//! It's based on the paper `Bit-parallel approximate string matching algorithms with transposition` from Heikki Hyyro
//!
//! ![benchmark results](https://raw.githubusercontent.com/maxbachmann/rapidfuzz-rs/main/rapidfuzz-benches/results/osa.svg)
//!
//! [`Damerau-Levenshtein`]: ../damerau_levenshtein/index.html

use crate::details::common::remove_common_affix;
use crate::details::distance::MetricUsize;
use crate::details::pattern_match_vector::{
    BitVectorInterface, BlockPatternMatchVector, PatternMatchVector,
};
use crate::HashableChar;
use std::marker::PhantomData;
use std::mem;

/// Bitparallel implementation of the OSA distance.
///
/// This implementation requires the first string to have a length <= 64.
/// The algorithm used stems from `hyrro_2002` and has a time complexity
/// of O(N). Comments and variable names in the implementation follow the
/// paper. This implementation is used internally when the strings are short enough
fn hyrroe2003<PmVec, Iter1, Iter2>(
    pm: &PmVec,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    _len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator,
    Iter2: Iterator,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    PmVec: BitVectorInterface,
{
    // VP is set to 1^m. Shifting by bitwidth would be undefined behavior
    let mut vp = !0_u64;
    let mut vn = 0_u64;
    let mut d0 = 0_u64;
    let mut pm_j_old = 0_u64;
    let mut curr_dist = len1;
    debug_assert!(len1 != 0);

    // mask used when computing D[m,j] in the paper 10^(m-1)
    let mask = 1_u64 << (len1 - 1);

    // Searching
    for ch2 in s2 {
        // Step 1: Computing D0
        let pm_j = pm.get(0, ch2);
        let tr = (((!d0) & pm_j) << 1) & pm_j_old;
        d0 = ((pm_j & vp).wrapping_add(vp) ^ vp) | pm_j | vn;
        d0 |= tr;

        // Step 2: Computing HP and HN
        let mut hp = vn | !(d0 | vp);
        let mut hn = d0 & vp;

        /* Step 3: Computing the value D[m,j] */
        curr_dist += usize::from(hp & mask != 0);
        curr_dist -= usize::from(hn & mask != 0);

        /* Step 4: Computing Vp and VN */
        hp = (hp << 1) | 1;
        hn <<= 1;

        vp = hn | !(d0 | hp);
        vn = hp & d0;
        pm_j_old = pm_j;
    }

    if curr_dist <= score_cutoff {
        Some(curr_dist)
    } else {
        None
    }
}

#[derive(Clone, Copy)]
struct OsaRow {
    vp: u64,
    vn: u64,
    d0: u64,
    pm: u64,
}

impl Default for OsaRow {
    fn default() -> Self {
        Self {
            vp: !0_u64,
            vn: 0_u64,
            d0: 0_u64,
            pm: 0_u64,
        }
    }
}

fn hyrroe2003_block<Iter1, Iter2>(
    pm: &BlockPatternMatchVector,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    _len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator,
    Iter2: Iterator,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    let word_size = 64;
    let words = pm.size();
    let last = 1_u64 << ((len1 - 1) % word_size);

    let mut curr_dist = len1;
    let mut old_vecs = vec![OsaRow::default(); words + 1];
    let mut new_vecs = vec![OsaRow::default(); words + 1];

    // Searching
    for ch2 in s2 {
        let mut hp_carry = 1_u64;
        let mut hn_carry = 0_u64;

        for word in 0..words {
            // retrieve bit vectors from last iterations
            let vn = old_vecs[word + 1].vn;
            let vp = old_vecs[word + 1].vp;
            let mut d0 = old_vecs[word + 1].d0;
            // D0 last word
            let d0_last = old_vecs[word].d0;

            // PM of last char same word
            let pm_j_old = old_vecs[word + 1].pm;
            // PM of last word
            let pm_last = new_vecs[word].pm;

            let pm_j = pm.get(word, ch2);
            let mut x = pm_j;
            let tr = ((((!d0) & x) << 1) | (((!d0_last) & pm_last) >> 63)) & pm_j_old;

            x |= hn_carry;
            d0 = ((x & vp).wrapping_add(vp) ^ vp) | x | vn | tr;

            let mut hp = vn | !(d0 | vp);
            let mut hn = d0 & vp;

            if word == words - 1 {
                curr_dist += usize::from(hp & last != 0);
                curr_dist -= usize::from(hn & last != 0);
            }

            let hp_carry_temp = hp_carry;
            hp_carry = hp >> 63;
            hp = (hp << 1) | hp_carry_temp;
            let hn_carry_temp = hn_carry;
            hn_carry = hn >> 63;
            hn = (hn << 1) | hn_carry_temp;

            new_vecs[word + 1].vp = hn | !(d0 | hp);
            new_vecs[word + 1].vn = hp & d0;
            new_vecs[word + 1].d0 = d0;
            new_vecs[word + 1].pm = pm_j;
        }

        mem::swap(&mut new_vecs, &mut old_vecs);
    }

    if curr_dist <= score_cutoff {
        Some(curr_dist)
    } else {
        None
    }
}

#[derive(Default, Clone)]
pub struct Dist;
#[derive(Default, Clone)]
pub struct Sim;
#[derive(Default, Clone)]
pub struct NormDist;
#[derive(Default, Clone)]
pub struct NormSim;
#[derive(Default, Clone)]
pub struct NoScoreCutoff;
#[derive(Default, Clone)]
pub struct WithScoreCutoff<T>(T);

pub struct IndividualComparator<T, ResultType, CutoffType> {
    score_cutoff: CutoffType,
    score_hint: Option<ResultType>,
    metric_type: PhantomData<T>,
}

impl<T, ResultType> Default for IndividualComparator<T, ResultType, NoScoreCutoff> {
    fn default() -> Self {
        Self {
            score_cutoff: Default::default(),
            score_hint: Default::default(),
            metric_type: PhantomData,
        }
    }
}

impl<T, ResultType, CutoffType> MetricUsize for IndividualComparator<T, ResultType, CutoffType> {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        if len1 < len2 {
            return self._distance(s2, len2, s1, len1, score_cutoff_, score_hint);
        }

        let score_cutoff = score_cutoff_.unwrap_or(usize::MAX);

        let affix = remove_common_affix(s1, len1, s2, len2);

        if affix.len1 == 0 {
            if affix.len2 <= score_cutoff {
                Some(affix.len2)
            } else {
                None
            }
        } else if affix.len1 <= 64 {
            let mut pm = PatternMatchVector::default();
            pm.insert(affix.s1.clone());
            hyrroe2003(
                &pm,
                affix.s1,
                affix.len1,
                affix.s2,
                affix.len2,
                score_cutoff,
            )
        } else {
            let mut pm = BlockPatternMatchVector::new(affix.len1);
            pm.insert(affix.s1.clone());
            hyrroe2003_block(
                &pm,
                affix.s1,
                affix.len1,
                affix.s2,
                affix.len2,
                score_cutoff,
            )
        }
    }
}

impl<T, ResultType, CutoffType> IndividualComparator<T, ResultType, CutoffType> {
    pub fn with_score_hint(mut self, score_hint: impl Into<ResultType>) -> Self {
        self.score_hint = Some(score_hint.into());
        self
    }

    pub fn with_score_cutoff(
        self,
        score_cutoff: impl Into<ResultType>,
    ) -> IndividualComparator<T, ResultType, WithScoreCutoff<ResultType>> {
        IndividualComparator {
            score_hint: self.score_hint,
            score_cutoff: WithScoreCutoff(score_cutoff.into()),
            metric_type: self.metric_type,
        }
    }
}

impl IndividualComparator<Dist, usize, NoScoreCutoff> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> usize
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.distance_(s1, s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl IndividualComparator<Sim, usize, NoScoreCutoff> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> usize
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.similarity_(s1, s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl IndividualComparator<NormDist, f64, NoScoreCutoff> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> f64
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.normalized_distance_(s1, s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl IndividualComparator<NormSim, f64, NoScoreCutoff> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> f64
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.normalized_similarity_(s1, s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl IndividualComparator<Dist, usize, WithScoreCutoff<usize>> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> Option<usize>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.distance_(s1, s2, Some(self.score_cutoff.0), self.score_hint)
    }
}

impl IndividualComparator<Sim, usize, WithScoreCutoff<usize>> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> Option<usize>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.similarity_(s1, s2, Some(self.score_cutoff.0), self.score_hint)
    }
}

impl IndividualComparator<NormDist, f64, WithScoreCutoff<f64>> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.normalized_distance_(s1, s2, Some(self.score_cutoff.0), self.score_hint)
    }
}

impl IndividualComparator<NormSim, f64, WithScoreCutoff<f64>> {
    pub fn compare<Iter1, Iter2>(&self, s1: Iter1, s2: Iter2) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        self.normalized_similarity_(s1, s2, Some(self.score_cutoff.0), self.score_hint)
    }
}

pub fn distance() -> IndividualComparator<Dist, usize, NoScoreCutoff> {
    IndividualComparator::<Dist, usize, NoScoreCutoff>::default()
}

pub fn similarity() -> IndividualComparator<Sim, usize, NoScoreCutoff> {
    IndividualComparator::<Sim, usize, NoScoreCutoff>::default()
}

pub fn normalized_distance() -> IndividualComparator<NormDist, f64, NoScoreCutoff> {
    IndividualComparator::<NormDist, f64, NoScoreCutoff>::default()
}

pub fn normalized_similarity() -> IndividualComparator<NormSim, f64, NoScoreCutoff> {
    IndividualComparator::<NormSim, f64, NoScoreCutoff>::default()
}

pub struct BatchComparator<Elem1> {
    s1: Vec<Elem1>,
    pm: BlockPatternMatchVector,
}

pub struct BatchComparatorBuilder<'a, Elem1, T, ResultType, CutoffType> {
    cache: &'a BatchComparator<Elem1>,
    score_cutoff: CutoffType,
    score_hint: Option<ResultType>,
    metric_type: PhantomData<T>,
}

impl<Elem1, T, ResultType, CutoffType> MetricUsize
    for BatchComparatorBuilder<'_, Elem1, T, ResultType, CutoffType>
{
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        _score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let score_cutoff = score_cutoff_.unwrap_or(usize::MAX);

        let dist = if self.cache.s1.is_empty() {
            len2
        } else if len2 == 0 {
            self.cache.s1.len()
        } else if self.cache.s1.len() <= 64 {
            hyrroe2003(&self.cache.pm, s1, len1, s2, len2, score_cutoff)?
        } else {
            hyrroe2003_block(&self.cache.pm, s1, len1, s2, len2, score_cutoff)?
        };

        if dist <= score_cutoff {
            Some(dist)
        } else {
            None
        }
    }
}

impl<'a, Elem1, T, ResultType> BatchComparatorBuilder<'a, Elem1, T, ResultType, NoScoreCutoff> {
    fn new(cache: &'a BatchComparator<Elem1>) -> Self {
        Self {
            cache,
            score_cutoff: Default::default(),
            score_hint: Default::default(),
            metric_type: PhantomData,
        }
    }
}

impl<'a, Elem1, T, ResultType, CutoffType>
    BatchComparatorBuilder<'a, Elem1, T, ResultType, CutoffType>
{
    pub fn with_score_hint(mut self, score_hint: impl Into<ResultType>) -> Self {
        self.score_hint = Some(score_hint.into());
        self
    }

    pub fn with_score_cutoff(
        self,
        score_cutoff: impl Into<ResultType>,
    ) -> BatchComparatorBuilder<'a, Elem1, T, ResultType, WithScoreCutoff<ResultType>> {
        BatchComparatorBuilder {
            cache: self.cache,
            score_hint: self.score_hint,
            score_cutoff: WithScoreCutoff(score_cutoff.into()),
            metric_type: self.metric_type,
        }
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, Dist, usize, NoScoreCutoff> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> usize
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        // todo this calls s1.count internally is this optimizes in the same way it
        // would if we pass in self.cache.s1.len()
        self.distance_(self.cache.s1.iter().copied(), s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, Sim, usize, NoScoreCutoff> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> usize
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.similarity_(self.cache.s1.iter().copied(), s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, NormDist, f64, NoScoreCutoff> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> f64
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.normalized_distance_(self.cache.s1.iter().copied(), s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, NormSim, f64, NoScoreCutoff> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> f64
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.normalized_similarity_(self.cache.s1.iter().copied(), s2, None, self.score_hint)
            .expect("no score_cutoff")
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, Dist, usize, WithScoreCutoff<usize>> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> Option<usize>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.distance_(
            self.cache.s1.iter().copied(),
            s2,
            Some(self.score_cutoff.0),
            self.score_hint,
        )
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, Sim, usize, WithScoreCutoff<usize>> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> Option<usize>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.similarity_(
            self.cache.s1.iter().copied(),
            s2,
            Some(self.score_cutoff.0),
            self.score_hint,
        )
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, NormDist, f64, WithScoreCutoff<f64>> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.normalized_distance_(
            self.cache.s1.iter().copied(),
            s2,
            Some(self.score_cutoff.0),
            self.score_hint,
        )
    }
}

impl<Elem1> BatchComparatorBuilder<'_, Elem1, NormSim, f64, WithScoreCutoff<f64>> {
    pub fn compare<Iter2>(&self, s2: Iter2) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.normalized_similarity_(
            self.cache.s1.iter().copied(),
            s2,
            Some(self.score_cutoff.0),
            self.score_hint,
        )
    }
}

impl<Elem1> BatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1_: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1_.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();

        let mut pm = BlockPatternMatchVector::new(s1.len());
        pm.insert(s1_iter);

        Self { s1, pm }
    }

    /// Normalized distance calculated similar to [`normalized_distance`]
    pub fn normalized_distance<'a>(
        &'a self,
    ) -> BatchComparatorBuilder<'a, Elem1, NormDist, f64, NoScoreCutoff> {
        BatchComparatorBuilder::<'a, Elem1, NormDist, f64, NoScoreCutoff>::new(self)
    }

    /// Normalized similarity calculated similar to [`normalized_similarity`]
    pub fn normalized_similarity<'a>(
        &'a self,
    ) -> BatchComparatorBuilder<'a, Elem1, NormSim, f64, NoScoreCutoff> {
        BatchComparatorBuilder::<'a, Elem1, NormSim, f64, NoScoreCutoff>::new(self)
    }

    /// Distance calculated similar to [`distance`]
    pub fn distance<'a>(&'a self) -> BatchComparatorBuilder<'a, Elem1, Dist, usize, NoScoreCutoff> {
        BatchComparatorBuilder::<'a, Elem1, Dist, usize, NoScoreCutoff>::new(self)
    }

    /// Similarity calculated similar to [`similarity`]
    pub fn similarity<'a>(
        &'a self,
    ) -> BatchComparatorBuilder<'a, Elem1, Sim, usize, NoScoreCutoff> {
        BatchComparatorBuilder::<'a, Elem1, Sim, usize, NoScoreCutoff>::new(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn _test_distance<Iter1, Iter2>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
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

        let scorer1 = distance()
            .with_score_hint(score_hint.unwrap_or(usize::MAX))
            .with_score_cutoff(score_cutoff.unwrap_or(usize::MAX));

        let res1 = scorer1.compare(s1.clone(), s2.clone());
        let res2 = scorer1.compare(s2.clone(), s1.clone());

        let scorer2 = BatchComparator::new(s1.clone());
        let res3 = scorer2
            .distance()
            .with_score_hint(score_hint.unwrap_or(usize::MAX))
            .with_score_cutoff(score_cutoff.unwrap_or(usize::MAX))
            .compare(s2.clone());
        let scorer3 = BatchComparator::new(s2.clone());
        let res4 = scorer3
            .distance()
            .with_score_hint(score_hint.unwrap_or(usize::MAX))
            .with_score_cutoff(score_cutoff.unwrap_or(usize::MAX))
            .compare(s1.clone());

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn _test_distance_ascii(
        s1: &str,
        s2: &str,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize> {
        let res1 = _test_distance(s1.chars(), s2.chars(), score_cutoff, score_hint);
        let res2 = _test_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    #[test]
    fn simple() {
        assert_eq!(Some(0), _test_distance_ascii("", "", None, None));

        assert_eq!(Some(4), _test_distance_ascii("aaaa", "", None, None));
        assert_eq!(None, _test_distance_ascii("aaaa", "", Some(1), None));

        assert_eq!(Some(3), _test_distance_ascii("CA", "ABC", None, None));
        assert_eq!(Some(1), _test_distance_ascii("CA", "AC", None, None));

        let filler = "a".repeat(64);
        let s1 = "a".to_string() + &filler + "CA" + &filler + "a";
        let s2 = "b".to_string() + &filler + "AC" + &filler + "b";
        assert_eq!(Some(3), _test_distance_ascii(&s1, &s2, None, None));
    }
}
