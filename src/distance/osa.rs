use crate::details::common::{
    find_common_prefix, find_common_suffix, norm_sim_to_norm_dist, HashableChar, UnrefIterator,
};
use crate::details::distance::{
    build_cached_distance_metric_funcs, build_cached_normalized_metric_funcs,
    build_distance_metric_funcs, build_normalized_metric_funcs,
};
use crate::details::pattern_match_vector::{
    BitVectorInterface, BitvectorHashmap, BlockPatternMatchVector, PatternMatchVector,
};
use std::mem;

/// Bitparallel implementation of the OSA distance.
///
/// This implementation requires the first string to have a length <= 64.
/// The algorithm used is described @cite hyrro_2002 and has a time complexity
/// of O(N). Comments and variable names in the implementation follow the
/// paper. This implementation is used internally when the strings are short enough
fn osa_hyrroe2003<PmVec, Iter1, Iter2, Elem1, Elem2>(
    pm: &PmVec,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    _len2: usize,
    score_cutoff: usize,
) -> usize
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
        curr_dist += (hp & mask != 0) as usize;
        curr_dist -= (hn & mask != 0) as usize;

        /* Step 4: Computing Vp and VN */
        hp = (hp << 1) | 1;
        hn <<= 1;

        vp = hn | !(d0 | hp);
        vn = hp & d0;
        pm_j_old = pm_j;
    }

    if curr_dist <= score_cutoff {
        curr_dist
    } else {
        score_cutoff + 1
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
        OsaRow {
            vp: !0_u64,
            vn: 0_u64,
            d0: 0_u64,
            pm: 0_u64,
        }
    }
}

fn osa_hyrroe2003_block<Iter1, Iter2, Elem1, Elem2>(
    pm: &BlockPatternMatchVector,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    _len2: usize,
    score_cutoff: usize,
) -> usize
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
    let mut old_vecs: Vec<OsaRow> = vec![Default::default(); words + 1];
    let mut new_vecs: Vec<OsaRow> = vec![Default::default(); words + 1];

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
                curr_dist += (hp & last != 0) as usize;
                curr_dist -= (hn & last != 0) as usize;
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
        curr_dist
    } else {
        score_cutoff + 1
    }
}

pub(crate) struct Osa {}

impl Osa {
    build_distance_metric_funcs!(Osa, usize, 0, usize::MAX);

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
        Elem1: PartialEq<Elem2>,
        Elem2: PartialEq<Elem1>,
    {
        len1.max(len2)
    }

    pub(crate) fn _distance<Iter1, Iter2, Elem1, Elem2>(
        s1: Iter1,
        mut len1: usize,
        s2: Iter2,
        mut len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
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
        if len1 < len2 {
            return Osa::_distance(s2, len2, s1, len1, score_cutoff, _score_hint);
        }

        let s1_iter_orig = s1.into_iter();
        let s2_iter_orig = s2.into_iter();
        let suffix_len = find_common_suffix(s1_iter_orig.clone(), s2_iter_orig.clone());
        let s1_iter_no_suffix = s1_iter_orig.take(len1 - suffix_len);
        let s2_iter_no_suffix = s2_iter_orig.take(len2 - suffix_len);
        let prefix_len = find_common_prefix(s1_iter_no_suffix.clone(), s2_iter_no_suffix.clone());
        let s1_iter = s1_iter_no_suffix.skip(prefix_len);
        let s2_iter = s2_iter_no_suffix.skip(prefix_len);
        len1 -= prefix_len + suffix_len;
        len2 -= prefix_len + suffix_len;

        if len1 == 0 {
            if len2 <= score_cutoff {
                len2
            } else {
                score_cutoff + 1
            }
        } else if len1 <= 64 {
            // rust fails to elide the copy when returning the array
            // from PatternMatchVector::new so manually inline it
            //let block = PatternMatchVector::new(s2_iter.clone());
            let mut pm = PatternMatchVector {
                map_unsigned: BitvectorHashmap::default(),
                map_signed: BitvectorHashmap::default(),
                extended_ascii: [0; 256],
            };
            pm.insert(s1_iter.clone());
            osa_hyrroe2003(&pm, s1_iter, len1, s2_iter, len2, score_cutoff)
        } else {
            let mut pm = BlockPatternMatchVector::new(len1);
            pm.insert(s1_iter.clone());
            osa_hyrroe2003_block(&pm, s1_iter, len1, s2_iter, len2, score_cutoff)
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
    Osa::distance(s1, s2, score_cutoff, score_hint)
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
    Osa::similarity(s1, s2, score_cutoff, score_hint)
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
    Osa::normalized_distance(s1, s2, score_cutoff, score_hint)
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
    Osa::normalized_similarity(s1, s2, score_cutoff, score_hint)
}

pub struct CachedOsa<Elem1>
where
    Elem1: HashableChar + Clone,
{
    s1: Vec<Elem1>,
    pm: BlockPatternMatchVector,
}

impl<Elem1> CachedOsa<Elem1>
where
    Elem1: HashableChar + Clone,
{
    build_cached_distance_metric_funcs!(CachedOsa, usize, 0, usize::MAX);

    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();

        let mut pm = BlockPatternMatchVector::new(s1.len());
        pm.insert(s1_iter);

        CachedOsa { s1, pm }
    }

    fn maximum<Iter2, Elem2>(&self, _s2: Iter2, len2: usize) -> usize
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
    {
        self.s1.len().max(len2)
    }

    fn _distance<Iter2, Elem2>(
        &self,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter2: IntoIterator<Item = Elem2>,
        Iter2::IntoIter: Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
        <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
    {
        let dist = if self.s1.is_empty() {
            len2
        } else if len2 == 0 {
            self.s1.len()
        } else if self.s1.len() <= 64 {
            osa_hyrroe2003(
                &self.pm,
                UnrefIterator {
                    seq: self.s1.iter(),
                },
                self.s1.len(),
                s2.into_iter(),
                len2,
                score_cutoff,
            )
        } else {
            osa_hyrroe2003_block(
                &self.pm,
                UnrefIterator {
                    seq: self.s1.iter(),
                },
                self.s1.len(),
                s2.into_iter(),
                len2,
                score_cutoff,
            )
        };

        if dist <= score_cutoff {
            dist
        } else {
            score_cutoff + 1
        }
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
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = distance(s1.clone(), s2.clone(), score_cutoff, score_hint);
        let res2 = distance(s2.clone(), s1.clone(), score_cutoff, score_hint);

        let scorer1 = CachedOsa::new(s1.clone());
        let res3 = scorer1.distance(s2.clone(), score_cutoff, score_hint);
        let scorer2 = CachedOsa::new(s2.clone());
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
    ) -> usize {
        let res1 = _test_distance(s1.chars(), s2.chars(), score_cutoff, score_hint);
        let res2 = _test_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    #[test]
    fn osa_simple() {
        assert_eq!(0, _test_distance_ascii("", "", None, None));

        assert_eq!(4, _test_distance_ascii("aaaa", "", None, None));
        assert_eq!(4, _test_distance_ascii("aaaa", "", None, None));
        assert_eq!(2, _test_distance_ascii("aaaa", "", Some(1), None));
        assert_eq!(2, _test_distance_ascii("aaaa", "", Some(1), None));

        assert_eq!(3, _test_distance_ascii("CA", "ABC", None, None));
        assert_eq!(1, _test_distance_ascii("CA", "AC", None, None));

        let filler = "a".repeat(64);
        let s1 = "a".to_string() + &filler + "CA" + &filler + "a";
        let s2 = "b".to_string() + &filler + "AC" + &filler + "b";
        assert_eq!(3, _test_distance_ascii(&s1, &s2, None, None));
    }
}
