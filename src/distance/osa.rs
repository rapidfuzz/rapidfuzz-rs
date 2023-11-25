use crate::details::common::{remove_common_affix, HashableChar};
use crate::details::distance::MetricUsize;
use crate::details::pattern_match_vector::{
    BitVectorInterface, BitvectorHashmap, BlockPatternMatchVector, PatternMatchVector,
};
use std::mem;

/// Bitparallel implementation of the OSA distance.
///
/// This implementation requires the first string to have a length <= 64.
/// The algorithm used is described @cite `hyrro_2002` and has a time complexity
/// of O(N). Comments and variable names in the implementation follow the
/// paper. This implementation is used internally when the strings are short enough
fn hyrroe2003<PmVec, Iter1, Iter2, Elem1, Elem2>(
    pm: &PmVec,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    _len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1>,
    Iter2: Iterator<Item = Elem2>,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
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

fn hyrroe2003_block<Iter1, Iter2, Elem1, Elem2>(
    pm: &BlockPatternMatchVector,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    _len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1>,
    Iter2: Iterator<Item = Elem2>,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
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

struct IndividualComparator;

impl MetricUsize for IndividualComparator {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
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
            // rust fails to elide the copy when returning the array
            // from PatternMatchVector::new so manually inline it
            //let block = PatternMatchVector::new(s2_iter.clone());
            let mut pm = PatternMatchVector {
                map_unsigned: BitvectorHashmap::default(),
                map_signed: BitvectorHashmap::default(),
                extended_ascii: [0; 256],
            };
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

pub struct BatchComparator<Elem1> {
    s1: Vec<Elem1>,
    pm: BlockPatternMatchVector,
}

impl<CharT> MetricUsize for BatchComparator<CharT> {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff_: Option<usize>,
        _score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let score_cutoff = score_cutoff_.unwrap_or(usize::MAX);

        let dist = if self.s1.is_empty() {
            len2
        } else if len2 == 0 {
            self.s1.len()
        } else if self.s1.len() <= 64 {
            hyrroe2003(&self.pm, s1, len1, s2, len2, score_cutoff)?
        } else {
            hyrroe2003_block(&self.pm, s1, len1, s2, len2, score_cutoff)?
        };

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
            self.s1.iter().copied(),
            self.s1.len(),
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
            self.s1.iter().copied(),
            self.s1.len(),
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
            self.s1.iter().copied(),
            self.s1.len(),
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

    fn _test_distance<Iter1, Iter2, Elem1, Elem2>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = distance(s1.clone(), s2.clone(), score_cutoff, score_hint);
        let res2 = distance(s2.clone(), s1.clone(), score_cutoff, score_hint);

        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.distance(s2.clone(), score_cutoff, score_hint);
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.distance(s1.clone(), score_cutoff, score_hint);

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
