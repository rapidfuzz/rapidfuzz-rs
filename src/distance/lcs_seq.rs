use crate::details::common::{remove_common_affix, HashableChar};
use crate::details::distance::MetricUsize;
use crate::details::intrinsics::{carrying_add, ceil_div_usize};
use crate::details::matrix::ShiftedBitMatrix;
use crate::details::pattern_match_vector::{
    BitVectorInterface, BlockPatternMatchVector, PatternMatchVector,
};
use std::cmp::{max, min};

#[derive(Default)]
struct ResultMatrix {
    s: ShiftedBitMatrix<u64>,
}

struct DistanceResult<const RECORD_MATRIX: usize> {
    record_matrix: [ResultMatrix; RECORD_MATRIX],
    sim: Option<usize>,
}

impl Default for DistanceResult<0> {
    fn default() -> Self {
        Self {
            record_matrix: [],
            sim: None,
        }
    }
}

impl Default for DistanceResult<1> {
    fn default() -> Self {
        Self {
            record_matrix: [ResultMatrix::default()],
            sim: None,
        }
    }
}

/// An encoded mbleven model table.
///
/// Each 8-bit integer represents an edit sequence, with using two
/// bits for a single operation.
///
/// Each Row of 8 integers represent all possible combinations
/// of edit sequences for a gived maximum edit distance and length
/// difference between the two strings, that is below the maximum
/// edit distance
///
///   0x1 = 01 = DELETE,
///   0x2 = 10 = INSERT
///
/// 0x5 -> DEL + DEL
/// 0x6 -> DEL + INS
/// 0x9 -> INS + DEL
/// 0xA -> INS + INS
static LCS_SEQ_MBLEVEN2018_MATRIX: [[u8; 6]; 14] = [
    // max edit distance 1
    // case does not occur
    [0x00, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 0
    [0x01, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 1
    // max edit distance 2
    [0x09, 0x06, 0x00, 0x00, 0x00, 0x00], // len_diff 0
    [0x01, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 1
    [0x05, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 2
    // max edit distance 3
    [0x09, 0x06, 0x00, 0x00, 0x00, 0x00], // len_diff 0
    [0x25, 0x19, 0x16, 0x00, 0x00, 0x00], // len_diff 1
    [0x05, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 2
    [0x15, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 3
    // max edit distance 4
    [0x96, 0x66, 0x5A, 0x99, 0x69, 0xA5], // len_diff 0
    [0x25, 0x19, 0x16, 0x00, 0x00, 0x00], // len_diff 1
    [0x65, 0x56, 0x95, 0x59, 0x00, 0x00], // len_diff 2
    [0x15, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 3
    [0x55, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 4
];

fn mbleven2018<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
{
    debug_assert!(len1 != 0);
    debug_assert!(len2 != 0);

    if len1 < len2 {
        return mbleven2018(s2, len2, s1, len1, score_cutoff);
    }

    let len_diff = len1 - len2;
    let max_misses = len1 + len2 - 2 * score_cutoff;
    let ops_index = (max_misses + max_misses * max_misses) / 2 + len_diff - 1;
    let possible_ops = &LCS_SEQ_MBLEVEN2018_MATRIX[ops_index];
    let mut max_len = 0;

    for &ops_ in possible_ops {
        let mut ops = ops_;
        let mut iter_s1 = s1.clone();
        let mut iter_s2 = s2.clone();
        let mut cur_len = 0;

        let mut cur1 = iter_s1.next();
        let mut cur2 = iter_s2.next();

        if ops == 0 {
            break;
        }

        while let (Some(ch1), Some(ch2)) = (&cur1, &cur2) {
            if ch1 == ch2 {
                cur_len += 1;
                cur1 = iter_s1.next();
                cur2 = iter_s2.next();
            } else {
                if ops == 0 {
                    break;
                }
                if (ops & 1) != 0 {
                    cur1 = iter_s1.next();
                } else if (ops & 2) != 0 {
                    cur2 = iter_s2.next();
                }

                ops >>= 2;
            }
        }

        max_len = max(max_len, cur_len);
    }

    if max_len >= score_cutoff {
        Some(max_len)
    } else {
        None
    }
}

fn lcs_unroll<const N: usize, const RECORD_MATRIX: usize, PmVec, Iter1, Iter2, Elem1, Elem2>(
    pm: &PmVec,
    _s1: Iter1,
    _len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> DistanceResult<RECORD_MATRIX>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    PmVec: BitVectorInterface,
    DistanceResult<RECORD_MATRIX>: Default,
{
    let mut s = [!0_u64; N];

    let mut res: DistanceResult<RECORD_MATRIX> = DistanceResult::default();
    if RECORD_MATRIX == 1 {
        res.record_matrix[0].s = ShiftedBitMatrix::<u64>::new(len2, N, !0_u64);
    }

    for (i, ch2) in s2.enumerate() {
        let mut carry = false;

        let mut calc = |word: usize| {
            let matches = pm.get(word, ch2);
            let u = s[word] & matches;
            let (x, carry_) = carrying_add(s[word], u, carry);
            carry = carry_;
            s[word] = x | (s[word] - u);

            if RECORD_MATRIX == 1 {
                *res.record_matrix[0].s.get_mut(i, word) = s[word];
            }
        };

        let unroll_factor = 3_usize;
        // manually unroll, since I am unaware of how to perform automatic loop unrolling in rust
        for j in 0..(N / unroll_factor) {
            calc(j * unroll_factor);
            calc((j + 1) * unroll_factor);
            calc((j + 2) * unroll_factor);
        }

        let j = N / unroll_factor;
        if N % unroll_factor == 2 {
            calc(j * unroll_factor);
            calc((j + 1) * unroll_factor);
        } else if N % unroll_factor == 1 {
            calc(j * unroll_factor);
        }
    }

    let mut sim = 0_usize;
    for x in s {
        sim += (!x).count_ones() as usize;
    }
    res.sim = if sim >= score_cutoff { Some(sim) } else { None };
    res
}

/// implementation is following the paper Bit-Parallel LCS-length Computation Revisited
/// from Heikki Hyyrö
///
/// The paper refers to s1 as m and s2 as n
fn lcs_blockwise<const RECORD_MATRIX: usize, PmVec, Iter1, Iter2, Elem1, Elem2>(
    pm: &PmVec,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> DistanceResult<RECORD_MATRIX>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    PmVec: BitVectorInterface,
    DistanceResult<RECORD_MATRIX>: Default,
{
    debug_assert!(score_cutoff <= len1);
    debug_assert!(score_cutoff <= len2);

    let word_size = 64;
    let words = pm.size();
    let mut s = vec![!0_u64; words];

    let band_width_left = len1 - score_cutoff;
    let band_width_right = len2 - score_cutoff;

    let mut res: DistanceResult<RECORD_MATRIX> = DistanceResult::default();
    if RECORD_MATRIX == 1 {
        let full_band = band_width_left + 1 + band_width_right;
        let full_band_words = min(words, full_band / word_size + 2);
        res.record_matrix[0].s = ShiftedBitMatrix::<u64>::new(len2, full_band_words, !0_u64);
    }

    // first_block is the index of the first block in Ukkonen band.
    let mut first_block = 0_usize;
    let mut last_block = min(words, ceil_div_usize(band_width_left + 1, word_size));

    for (row, ch2) in s2.enumerate() {
        let mut carry = false;

        if RECORD_MATRIX == 1 {
            res.record_matrix[0]
                .s
                .set_offset(row, (first_block * word_size) as isize);
        }

        for (word, s_ref) in s.iter_mut().enumerate().take(last_block).skip(first_block) {
            let matches = pm.get(word, ch2);
            let u = *s_ref & matches;
            let (x, carry_) = carrying_add(*s_ref, u, carry);
            carry = carry_;
            *s_ref = x | (*s_ref - u);

            if RECORD_MATRIX == 1 {
                *res.record_matrix[0].s.get_mut(row, word - first_block) = *s_ref;
            }
        }

        if row > band_width_right {
            first_block = (row - band_width_right) / word_size;
        }

        if row + 1 + band_width_left <= len1 {
            last_block = ceil_div_usize(row + 1 + band_width_left, word_size);
        }
    }

    let mut sim = 0_usize;
    for x in s {
        sim += (!x).count_ones() as usize;
    }
    res.sim = if sim >= score_cutoff { Some(sim) } else { None };
    res
}

fn longest_common_subsequence_with_pm<PmVec, Iter1, Iter2, Elem1, Elem2>(
    pm: &PmVec,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    PmVec: BitVectorInterface,
{
    let word_size = 64;
    let words = pm.size();
    let band_width_left = len1 - score_cutoff;
    let band_width_right = len2 - score_cutoff;
    let full_band = band_width_left + 1 + band_width_right;
    let full_band_words = min(words, (full_band / word_size) + 2);

    if full_band_words < words {
        let res: DistanceResult<0> = lcs_blockwise(pm, s1, len1, s2, len2, score_cutoff);
        return res.sim;
    }

    match ceil_div_usize(len1, word_size) {
        0 => {
            if 0 == score_cutoff {
                Some(0)
            } else {
                None
            }
        }
        1 => {
            let func = lcs_unroll::<1, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        2 => {
            let func = lcs_unroll::<2, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        3 => {
            let func = lcs_unroll::<3, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        4 => {
            let func = lcs_unroll::<4, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        5 => {
            let func = lcs_unroll::<5, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        6 => {
            let func = lcs_unroll::<6, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        7 => {
            let func = lcs_unroll::<7, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        8 => {
            let func = lcs_unroll::<8, 0, PmVec, Iter1, Iter2, Elem1, Elem2>;
            func(pm, s1, len1, s2, len2, score_cutoff).sim
        }
        _ => {
            let res: DistanceResult<0> = lcs_blockwise(pm, s1, len1, s2, len2, score_cutoff);
            res.sim
        }
    }
}

fn longest_common_subsequence_without_pm<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    if len1 == 0 {
        Some(0)
    } else if len1 <= 64 {
        let mut pm = PatternMatchVector::default();
        pm.insert(s1.clone());
        longest_common_subsequence_with_pm(&pm, s1, len1, s2, len2, score_cutoff)
    } else {
        // todo add score_hint support
        let mut pm = BlockPatternMatchVector::new(len1);
        pm.insert(s1.clone());

        longest_common_subsequence_with_pm(&pm, s1, len1, s2, len2, score_cutoff)
    }
}

pub(crate) fn similarity_with_pm<PmVec, Iter1, Iter2, Elem1, Elem2>(
    pm: &PmVec,
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
    PmVec: BitVectorInterface,
{
    if score_cutoff > len1 || score_cutoff > len2 {
        return None;
    }

    let max_misses = len1 + len2 - 2 * score_cutoff;
    // no edits are allowed
    if max_misses == 0 || (max_misses == 1 && len1 == len2) {
        return if s1.into_iter().eq(s2) {
            Some(len1)
        } else {
            None
        };
    }

    if max_misses < len1.abs_diff(len2) {
        return None;
    }

    // do this first, since we can not remove any affix in encoded form
    if max_misses >= 5 {
        return longest_common_subsequence_with_pm(pm, s1, len1, s2, len2, score_cutoff);
    }

    // remove common affix and count it as part of the LCS
    let affix = remove_common_affix(s1, len1, s2, len2);
    let mut lcs_sim = affix.prefix_len + affix.suffix_len;
    if len1 != 0 && len2 != 0 {
        lcs_sim += mbleven2018(
            affix.s1,
            affix.len1,
            affix.s2,
            affix.len2,
            score_cutoff - lcs_sim,
        )?;
    }

    if lcs_sim >= score_cutoff {
        Some(lcs_sim)
    } else {
        None
    }
}

fn similarity_without_pm<Iter1, Iter2, Elem1, Elem2>(
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
    // Swapping the strings so the second string is shorter
    if len1 < len2 {
        return similarity_without_pm(s2, len2, s1, len1, score_cutoff);
    }

    if score_cutoff > len1 || score_cutoff > len2 {
        return None;
    }

    let max_misses = len1 + len2 - 2 * score_cutoff;
    // no edits are allowed
    if max_misses == 0 || (max_misses == 1 && len1 == len2) {
        return if s1.into_iter().eq(s2) {
            Some(len1)
        } else {
            None
        };
    }

    if max_misses < len1.abs_diff(len2) {
        return None;
    }

    // remove common affix and count it as part of the LCS
    let affix = remove_common_affix(s1, len1, s2, len2);

    let mut lcs_sim = affix.prefix_len + affix.suffix_len;
    if len1 != 0 && len2 != 0 {
        let adjusted_cutoff = if score_cutoff >= lcs_sim {
            score_cutoff - lcs_sim
        } else {
            0
        };
        if max_misses < 5 {
            lcs_sim += mbleven2018(affix.s1, affix.len1, affix.s2, affix.len2, adjusted_cutoff)?;
        } else {
            lcs_sim += longest_common_subsequence_without_pm(
                affix.s1,
                affix.len1,
                affix.s2,
                affix.len2,
                adjusted_cutoff,
            )?;
        }
    }

    if lcs_sim >= score_cutoff {
        Some(lcs_sim)
    } else {
        None
    }
}

pub(crate) struct IndividualComparator;

impl MetricUsize for IndividualComparator {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        max(len1, len2)
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        _score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        similarity_without_pm(s1, len1, s2, len2, score_cutoff.unwrap_or(0))
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
    pub(crate) s1: Vec<Elem1>,
    pub(crate) pm: BlockPatternMatchVector,
}

impl<CharT> MetricUsize for BatchComparator<CharT> {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        len1.max(len2)
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        _score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        similarity_with_pm(&self.pm, s1, len1, s2, len2, score_cutoff.unwrap_or(0))
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

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            match ($x, $y) {
                (None, None) => {}
                (Some(val1), Some(val2)) => {
                    if (val1 - val2).abs() > $d {
                        panic!("{:?} != {:?}", $x, $y);
                    }
                }
                (_, _) => panic!("{:?} != {:?}", $x, $y),
            }
        };
    }

    fn test_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
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
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = distance(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = distance(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );

        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.distance(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.distance(s1.clone(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn test_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
    where
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let res1 = test_distance(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    fn test_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
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
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = similarity(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = similarity(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );

        let scorer1 = BatchComparator::new(s1.clone());
        let res3 = scorer1.similarity(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.similarity(s1.clone(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn test_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
    where
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let res1 = test_similarity(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    fn test_normalized_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
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
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = normalized_distance(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = normalized_distance(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let scorer1 = BatchComparator::new(s1.clone());
        let res3 =
            scorer1.normalized_distance(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.normalized_distance(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn test_normalized_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = test_normalized_distance(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_normalized_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    fn test_normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
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
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let s1 = s1_.into_iter();
        let s2 = s2_.into_iter();
        let res1 = normalized_similarity(
            s1.clone(),
            s2.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = normalized_similarity(
            s2.clone(),
            s1.clone(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let scorer1 = BatchComparator::new(s1.clone());
        let res3 =
            scorer1.normalized_similarity(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone());
        let res4 = scorer2.normalized_similarity(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn test_normalized_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = test_normalized_similarity(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = test_normalized_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    #[test]
    fn similar() {
        assert_eq!(Some(0), test_distance_ascii("aaaa", "aaaa", None, None));
        assert_eq!(Some(4), test_similarity_ascii("aaaa", "aaaa", None, None));
        assert_delta!(
            Some(0.0),
            test_normalized_distance_ascii("aaaa", "aaaa", None, None),
            0.0001
        );
        assert_delta!(
            Some(1.0),
            test_normalized_similarity_ascii("aaaa", "aaaa", None, None),
            0.0001
        );
    }

    #[test]
    fn completely_different() {
        assert_eq!(Some(4), test_distance_ascii("aaaa", "bbbb", None, None));
        assert_eq!(Some(0), test_similarity_ascii("aaaa", "bbbb", None, None));
        assert_delta!(
            Some(1.0),
            test_normalized_distance_ascii("aaaa", "bbbb", None, None),
            0.0001
        );
        assert_delta!(
            Some(0.0),
            test_normalized_similarity_ascii("aaaa", "bbbb", None, None),
            0.0001
        );
    }

    /// test mbleven implementation
    #[test]
    fn test_mbleven() {
        let mut a = "South Korea";
        let mut b = "North Korea";

        assert_eq!(Some(9), test_similarity_ascii(a, b, None, None));
        assert_eq!(Some(9), test_similarity_ascii(a, b, 9, None));
        assert_eq!(None, test_similarity_ascii(a, b, 10, None));

        assert_eq!(Some(2), test_distance_ascii(a, b, None, None));
        assert_eq!(Some(2), test_distance_ascii(a, b, 4, None));
        assert_eq!(Some(2), test_distance_ascii(a, b, 3, None));
        assert_eq!(Some(2), test_distance_ascii(a, b, 2, None));
        assert_eq!(None, test_distance_ascii(a, b, 1, None));
        assert_eq!(None, test_distance_ascii(a, b, 0, None));

        a = "aabc";
        b = "cccd";
        assert_eq!(Some(1), test_similarity_ascii(a, b, None, None));
        assert_eq!(Some(1), test_similarity_ascii(a, b, 1, None));
        assert_eq!(None, test_similarity_ascii(a, b, 2, None));

        assert_eq!(Some(3), test_distance_ascii(a, b, None, None));
        assert_eq!(Some(3), test_distance_ascii(a, b, 4, None));
        assert_eq!(Some(3), test_distance_ascii(a, b, 3, None));
        assert_eq!(None, test_distance_ascii(a, b, 2, None));
        assert_eq!(None, test_distance_ascii(a, b, 1, None));
        assert_eq!(None, test_distance_ascii(a, b, 0, None));
    }

    #[test]
    fn test_cached() {
        let a = "001";
        let b = "220";
        assert_eq!(Some(1), test_similarity_ascii(a, b, None, None));
    }
}
