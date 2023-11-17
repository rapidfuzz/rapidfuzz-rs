use crate::details::common::{
    find_common_prefix, find_common_suffix, norm_sim_to_norm_dist, HashableChar, UnrefIterator,
};
use crate::details::distance::{
    build_cached_distance_metric_funcs, build_cached_normalized_metric_funcs,
    build_distance_metric_funcs, build_normalized_metric_funcs,
};
use crate::details::growing_hashmap::HybridGrowingHashmap;
use std::cmp::{max, min};
use std::mem;

#[derive(Clone, Copy, PartialEq, Eq)]
struct RowId {
    val: isize,
}

impl Default for RowId {
    fn default() -> Self {
        RowId { val: -1 }
    }
}

/// based on the paper
/// "Linear space string correction algorithm using the Damerau-Levenshtein distance"
/// from Chunchun Zhao and Sartaj Sahni
///
/// todo in c++ this is templated on an integer type which reduced
/// memory usage depending on string lengths
fn damerau_damerau_levenshtein_distance_zhao<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> usize
where
    Iter1: Iterator<Item = Elem1>,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
{
    let max_val = max(len1, len2) as isize + 1;

    let mut last_row_id = HybridGrowingHashmap::<RowId> {
        map_unsigned: Default::default(),
        map_signed: Default::default(),
        extended_ascii: [Default::default(); 256],
    };
    let size = len2 + 2;
    let mut fr = vec![max_val; size];
    let mut r1 = vec![max_val; size];
    let mut r: Vec<isize> = (max_val..max_val + 1)
        .chain(0..(size - 1) as isize)
        .collect();

    for (i, ch1) in s1.enumerate().map(|(i, ch1)| (i + 1, ch1)) {
        mem::swap(&mut r, &mut r1);
        let mut last_col_id: isize = -1;
        let mut last_i2l1 = r[1];
        r[1] = i as isize;
        let mut t = max_val;

        for (j, ch2) in s2.clone().enumerate().map(|(j, ch2)| (j + 1, ch2)) {
            let diag = r1[j] + (ch1 != ch2) as isize;
            let left = r[j] + 1;
            let up = r1[j + 1] + 1;
            let mut temp = min(diag, min(left, up));

            if ch1 == ch2 {
                last_col_id = j as isize; // last occurence of s1_i
                fr[j + 1] = r1[j - 1]; // save H_k-1,j-2
                t = last_i2l1; // save H_i-2,l-1
            } else {
                let k = last_row_id.get(ch2).val;
                let l = last_col_id;

                if j as isize - l == 1 {
                    let transpose = fr[j + 1] + (i as isize - k);
                    temp = min(temp, transpose);
                } else if i as isize - k == 1 {
                    let transpose = t + (j as isize - l);
                    temp = min(temp, transpose);
                }
            }

            last_i2l1 = r[j + 1];
            r[j + 1] = temp;
        }

        last_row_id.get_mut(ch1).val = i as isize;
    }

    let dist = r[len2 + 1] as usize;

    if dist < score_cutoff {
        dist
    } else {
        score_cutoff + 1
    }
}

fn damerau_damerau_levenshtein_distance_impl<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    mut len1: usize,
    s2: Iter2,
    mut len2: usize,
    score_cutoff: usize,
) -> usize
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
{
    if score_cutoff < len1.abs_diff(len2) {
        return score_cutoff + 1;
    }

    // common affix does not effect Levenshtein distance
    let suffix_len = find_common_suffix(s1.clone(), s2.clone());
    let s1_iter_no_suffix = s1.take(len1 - suffix_len);
    let s2_iter_no_suffix = s2.take(len2 - suffix_len);
    let prefix_len = find_common_prefix(s1_iter_no_suffix.clone(), s2_iter_no_suffix.clone());
    let s1_iter = s1_iter_no_suffix.skip(prefix_len);
    let s2_iter = s2_iter_no_suffix.skip(prefix_len);
    len1 -= prefix_len + suffix_len;
    len2 -= prefix_len + suffix_len;

    damerau_damerau_levenshtein_distance_zhao(s1_iter, len1, s2_iter, len2, score_cutoff)
}

pub(crate) struct DamerauLevenshtein {}

impl DamerauLevenshtein {
    build_distance_metric_funcs!(DamerauLevenshtein, usize, 0, usize::MAX);

    fn maximum(len1: usize, len2: usize) -> usize {
        max(len1, len2)
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        damerau_damerau_levenshtein_distance_impl(s1, len1, s2, len2, score_cutoff)
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
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    DamerauLevenshtein::distance(s1, s2, score_cutoff, score_hint)
}

pub fn similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: Option<usize>,
    score_hint: Option<usize>,
) -> usize
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    DamerauLevenshtein::similarity(s1, s2, score_cutoff, score_hint)
}

pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: Option<f64>,
    score_hint: Option<f64>,
) -> f64
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    DamerauLevenshtein::normalized_distance(s1, s2, score_cutoff, score_hint)
}

pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: Option<f64>,
    score_hint: Option<f64>,
) -> f64
where
    Iter1: IntoIterator<Item = Elem1>,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator<Item = Elem2>,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    DamerauLevenshtein::normalized_similarity(s1, s2, score_cutoff, score_hint)
}

pub struct CachedDamerauLevenshtein<Elem1>
where
    Elem1: HashableChar + Clone,
{
    s1: Vec<Elem1>,
}

impl<Elem1> CachedDamerauLevenshtein<Elem1>
where
    Elem1: HashableChar + Clone,
{
    build_cached_distance_metric_funcs!(CachedDamerauLevenshtein, usize, 0, usize::MAX);

    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();
        CachedDamerauLevenshtein { s1 }
    }

    fn maximum(&self, len2: usize) -> usize {
        max(self.s1.len(), len2)
    }

    fn _distance<Iter2, Elem2>(
        &self,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        _score_hint: usize,
    ) -> usize
    where
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        damerau_damerau_levenshtein_distance_impl(
            UnrefIterator {
                seq: self.s1.iter(),
            },
            self.s1.len(),
            s2,
            len2,
            score_cutoff,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static EMPTY: &str = "";
    static TEST: &str = "aaaa";
    static NO_SUFFIX: &str = "aaa";
    static NO_SUFFIX2: &str = "aaab";
    static SWAPPED1: &str = "abaa";
    static SWAPPED2: &str = "baaa";
    static REPLACE_ALL: &str = "bbbb";

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            if ($x - $y).abs() > $d {
                panic!();
            }
        };
    }

    fn _test_distance<Iter1, Iter2, Elem1, Elem2>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> usize
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

        let scorer1 = CachedDamerauLevenshtein::new(s1.clone());
        let res3 = scorer1.distance(s2.clone(), score_cutoff, score_hint);
        let scorer2 = CachedDamerauLevenshtein::new(s2.clone());
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

    fn _test_normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> f64
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
        let res1 = normalized_similarity(s1.clone(), s2.clone(), score_cutoff, score_hint);
        let res2 = normalized_similarity(s2.clone(), s1.clone(), score_cutoff, score_hint);
        let scorer1 = CachedDamerauLevenshtein::new(s1.clone());
        let res3 = scorer1.normalized_similarity(s2.clone(), score_cutoff, score_hint);
        let scorer2 = CachedDamerauLevenshtein::new(s2.clone());
        let res4 = scorer2.normalized_similarity(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn _test_normalized_similarity_ascii(
        s1: &str,
        s2: &str,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> f64 {
        let res1 = _test_normalized_similarity(s1.chars(), s2.chars(), score_cutoff, score_hint);
        let res2 = _test_normalized_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    /// levenshtein calculates empty sequence
    #[test]
    fn damerau_levenshtein_empty() {
        assert_eq!(0, _test_distance_ascii(EMPTY, EMPTY, None, None));
        assert_eq!(4, _test_distance_ascii(TEST, EMPTY, None, None));
        assert_eq!(4, _test_distance_ascii(EMPTY, TEST, None, None));
    }

    /// levenshtein calculates correct distances
    #[test]
    fn damerau_levenshtein_simple() {
        assert_eq!(0, _test_distance_ascii(TEST, TEST, None, None));
        assert_eq!(1, _test_distance_ascii(TEST, NO_SUFFIX, None, None));
        assert_eq!(1, _test_distance_ascii(TEST, NO_SUFFIX2, None, None));
        assert_eq!(1, _test_distance_ascii(SWAPPED1, SWAPPED2, None, None));
        assert_eq!(4, _test_distance_ascii(TEST, REPLACE_ALL, None, None));

        assert_delta!(
            1.0,
            _test_normalized_similarity_ascii(TEST, TEST, None, None),
            0.0001
        );
        assert_delta!(
            0.75,
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX, None, None),
            0.0001
        );
        assert_delta!(
            0.75,
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX2, None, None),
            0.0001
        );
        assert_delta!(
            0.75,
            _test_normalized_similarity_ascii(SWAPPED1, SWAPPED2, None, None),
            0.0001
        );
        assert_delta!(
            0.0,
            _test_normalized_similarity_ascii(TEST, REPLACE_ALL, None, None),
            0.0001
        );
    }
}
