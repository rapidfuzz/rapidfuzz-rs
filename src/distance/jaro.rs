//! Jaro similarity
//!
//! The Jaro similarity is a measure of similarity between two strings, often
//! used in the field of record linkage and string matching. It's particularly
//! effective in comparing short strings, such as names.
//! The algorithm considers both the common characters and their order in the strings,
//! as well as the number of transpositions needed to make the strings equal.
//!
//! # Performance
//!
//! The implementation has a runtime complexity of `O([N/64]*M)` and a memory usage of `O(N)`.
//!
//! ![benchmark results](https://raw.githubusercontent.com/maxbachmann/rapidfuzz-rs/main/rapidfuzz-benches/results/jaro.svg)
//!

use crate::details::common::find_common_prefix;
use crate::details::distance::Metricf64;
use crate::details::intrinsics::{bit_mask_lsb_u64, blsi_u64, ceil_div_usize};
use crate::details::pattern_match_vector::{
    BitVectorInterface, BlockPatternMatchVector, PatternMatchVector,
};
use crate::Hash;
use crate::HashableChar;
use std::cmp::min;

struct FlaggedCharsWord {
    p_flag: u64,
    t_flag: u64,
}

impl FlaggedCharsWord {
    const fn count_common_chars(&self) -> usize {
        self.p_flag.count_ones() as usize
    }
}

struct FlaggedCharsMultiword {
    p_flag: Vec<u64>,
    t_flag: Vec<u64>,
}

impl FlaggedCharsMultiword {
    fn count_common_chars(&self) -> usize {
        if self.p_flag.len() < self.t_flag.len() {
            self.p_flag.iter().map(|x| x.count_ones() as usize).sum()
        } else {
            self.t_flag.iter().map(|x| x.count_ones() as usize).sum()
        }
    }
}

struct SearchBoundMask {
    words: usize,
    empty_words: usize,
    last_mask: u64,
    first_mask: u64,
}

fn calculate_similarity(
    p_len: usize,
    t_len: usize,
    common_chars: usize,
    mut transposition: usize,
) -> f64 {
    transposition /= 2;
    let mut sim: f64 = 0.0;
    sim += common_chars as f64 / p_len as f64;
    sim += common_chars as f64 / t_len as f64;
    sim += (common_chars as f64 - transposition as f64) / common_chars as f64;

    sim / 3.0
}

// filter matches below score_cutoff based on string lengths
fn length_filter(p_len: usize, t_len: usize, score_cutoff: f64) -> bool {
    if t_len == 0 || p_len == 0 {
        return false;
    }

    let min_len = p_len.min(t_len) as f64;
    let mut sim = min_len / p_len as f64 + min_len / t_len as f64 + 1.0;
    sim /= 3.0;
    sim >= score_cutoff
}

// filter matches below score_cutoff based on string lengths and common characters
fn common_char_filter(p_len: usize, t_len: usize, common_chars: usize, score_cutoff: f64) -> bool {
    if common_chars == 0 {
        return false;
    }

    let mut sim: f64 = 0.0;
    sim += common_chars as f64 / p_len as f64;
    sim += common_chars as f64 / t_len as f64;
    sim += 1.0;
    sim /= 3.0;
    sim >= score_cutoff
}

fn flag_similar_characters_word<PmVec, Iter2>(
    pm: &PmVec,
    len1: usize,
    mut s2: Iter2,
    len2: usize,
    bound: usize,
) -> FlaggedCharsWord
where
    Iter2: Iterator,
    Iter2::Item: HashableChar + Copy,
    PmVec: BitVectorInterface,
{
    debug_assert!(len1 <= 64);
    debug_assert!(len2 <= 64);
    debug_assert!(bound > len1 || len1 - bound <= len2);

    let mut flagged = FlaggedCharsWord {
        p_flag: 0,
        t_flag: 0,
    };

    let mut bound_mask = bit_mask_lsb_u64(bound + 1);

    let mut j = 0;
    for ch2 in (&mut s2).take(bound) {
        let pm_j = pm.get(0, ch2) & bound_mask & !flagged.p_flag;
        flagged.p_flag |= blsi_u64(pm_j);
        flagged.t_flag |= u64::from(pm_j != 0) << j;

        bound_mask = (bound_mask << 1) | 1;
        j += 1;
    }

    for ch2 in s2 {
        let pm_j = pm.get(0, ch2) & bound_mask & !flagged.p_flag;
        flagged.p_flag |= blsi_u64(pm_j);
        flagged.t_flag |= u64::from(pm_j != 0) << j;

        bound_mask <<= 1;
        j += 1;
    }

    flagged
}

fn flag_similar_characters_step<CharT>(
    pm: &BlockPatternMatchVector,
    t_j: CharT,
    flagged: &mut FlaggedCharsMultiword,
    j: usize,
    bound_mask: &SearchBoundMask,
) where
    CharT: HashableChar + Copy,
{
    let j_word = j / 64;
    let j_pos = j % 64;
    let mut word = bound_mask.empty_words;
    let last_word = word + bound_mask.words;

    if bound_mask.words == 1 {
        let pm_j = pm.get(word, t_j)
            & bound_mask.last_mask
            & bound_mask.first_mask
            & (!flagged.p_flag[word]);

        flagged.p_flag[word] |= blsi_u64(pm_j);
        flagged.t_flag[j_word] |= u64::from(pm_j != 0) << j_pos;
        return;
    }

    if bound_mask.first_mask != 0 {
        let pm_j = pm.get(word, t_j) & bound_mask.first_mask & (!flagged.p_flag[word]);

        if pm_j != 0 {
            flagged.p_flag[word] |= blsi_u64(pm_j);
            flagged.t_flag[j_word] |= 1_u64 << j_pos;
            return;
        }
        word += 1;
    }

    // unroll for better performance on long sequences when access is fast
    let is_ascii = match t_j.hash_char() {
        Hash::UNSIGNED(value) => value < 256,
        Hash::SIGNED(value) => value < 256,
    };
    if is_ascii {
        while word + 3 < last_word - 1 {
            let pm_j = [
                pm.get(word, t_j) & (!flagged.p_flag[word]),
                pm.get(word + 1, t_j) & (!flagged.p_flag[word + 1]),
                pm.get(word + 2, t_j) & (!flagged.p_flag[word + 2]),
                pm.get(word + 3, t_j) & (!flagged.p_flag[word + 3]),
            ];

            if pm_j[0] != 0 {
                flagged.p_flag[word] |= blsi_u64(pm_j[0]);
                flagged.t_flag[j_word] |= 1_u64 << j_pos;
                return;
            }
            if pm_j[1] != 0 {
                flagged.p_flag[word] |= blsi_u64(pm_j[1]);
                flagged.t_flag[j_word] |= 1_u64 << j_pos;
                return;
            }
            if pm_j[2] != 0 {
                flagged.p_flag[word] |= blsi_u64(pm_j[2]);
                flagged.t_flag[j_word] |= 1_u64 << j_pos;
                return;
            }
            if pm_j[3] != 0 {
                flagged.p_flag[word] |= blsi_u64(pm_j[3]);
                flagged.t_flag[j_word] |= 1_u64 << j_pos;
                return;
            }

            word += 3;
        }
    }

    while word < last_word - 1 {
        let pm_j = pm.get(word, t_j) & (!flagged.p_flag[word]);

        if pm_j != 0 {
            flagged.p_flag[word] |= blsi_u64(pm_j);
            flagged.t_flag[j_word] |= 1_u64 << j_pos;
            return;
        }
        word += 1;
    }

    if bound_mask.last_mask != 0 {
        let pm_j = pm.get(word, t_j) & bound_mask.last_mask & (!flagged.p_flag[word]);

        flagged.p_flag[word] |= blsi_u64(pm_j);
        flagged.t_flag[j_word] |= u64::from(pm_j != 0) << j_pos;
    }
}

fn flag_similar_characters_block<Iter2>(
    pm: &BlockPatternMatchVector,
    len1: usize,
    s2: Iter2,
    len2: usize,
    bound: usize,
) -> FlaggedCharsMultiword
where
    Iter2: Iterator,
    Iter2::Item: HashableChar + Copy,
{
    debug_assert!(len1 > 64 || len2 > 64);
    debug_assert!(bound > len1 || len1 - bound <= len2);
    debug_assert!(bound >= 31);

    let mut flagged = FlaggedCharsMultiword {
        p_flag: vec![0; ceil_div_usize(len1, 64)],
        t_flag: vec![0; ceil_div_usize(len2, 64)],
    };

    let start_range = min(bound + 1, len1);
    let mut bound_mask = SearchBoundMask {
        words: 1 + start_range / 64,
        empty_words: 0,
        last_mask: (1_u64 << (start_range % 64)) - 1,
        first_mask: !0_u64,
    };

    for (j, ch2) in s2.enumerate() {
        flag_similar_characters_step(pm, ch2, &mut flagged, j, &bound_mask);

        if j + bound + 1 < len1 {
            bound_mask.last_mask = (bound_mask.last_mask << 1) | 1;
            if j + bound + 2 < len1 && bound_mask.last_mask == !0_u64 {
                bound_mask.last_mask = 0;
                bound_mask.words += 1;
            }
        }

        if j >= bound {
            bound_mask.first_mask <<= 1;
            if bound_mask.first_mask == 0 {
                bound_mask.first_mask = !0_u64;
                bound_mask.words -= 1;
                bound_mask.empty_words += 1;
            }
        }
    }

    flagged
}

fn count_transpositions_word<PmVec, Iter2>(
    pm: &PmVec,
    mut s2: Iter2,
    _len2: usize,
    flagged: &FlaggedCharsWord,
) -> usize
where
    Iter2: Iterator + Clone,
    PmVec: BitVectorInterface,
    Iter2::Item: HashableChar,
{
    let mut p_flag = flagged.p_flag;
    let mut t_flag = flagged.t_flag;
    let mut transpositions = 0_usize;
    while t_flag != 0 {
        let pattern_flag_mask = blsi_u64(p_flag);

        let s2_index = t_flag.trailing_zeros() as usize;
        let ch2 = s2
            .nth(s2_index)
            .expect("these can't be outside, since we set the flags based on available indexes");

        transpositions += usize::from((pm.get(0, ch2) & pattern_flag_mask) == 0);

        t_flag = (t_flag >> 1) >> s2_index;
        p_flag ^= pattern_flag_mask;
    }

    transpositions
}

fn count_transpositions_block<Iter2>(
    pm: &BlockPatternMatchVector,
    mut s2: Iter2,
    _len2: usize,
    flagged: &FlaggedCharsMultiword,
    mut flagged_chars: usize,
) -> usize
where
    Iter2: Iterator + Clone,
    Iter2::Item: HashableChar,
{
    let mut text_word: usize = 0;
    let mut pattern_word: usize = 0;
    let mut t_flag = flagged.t_flag[text_word];
    let mut p_flag = flagged.p_flag[pattern_word];

    let mut transpositions = 0;
    let mut s2_pos = 0_usize;
    while flagged_chars != 0 {
        while t_flag == 0 {
            text_word += 1;
            if s2_pos < 64 {
                s2.nth(64 - 1 - s2_pos);
            }
            t_flag = flagged.t_flag[text_word];
            s2_pos = 0;
        }
        while t_flag != 0 {
            while p_flag == 0 {
                pattern_word += 1;
                p_flag = flagged.p_flag[pattern_word];
            }

            let pattern_flag_mask = blsi_u64(p_flag);
            let s2_index = t_flag.trailing_zeros() as usize;
            let ch2 = s2.nth(s2_index).expect(
                "these can't be outside, since we set the flags based on available indexes",
            );
            s2_pos += s2_index + 1;

            transpositions += usize::from((pm.get(pattern_word, ch2) & pattern_flag_mask) == 0);

            t_flag = (t_flag >> 1) >> s2_index;
            p_flag ^= pattern_flag_mask;
            flagged_chars -= 1;
        }
    }

    transpositions
}

pub(crate) fn similarity_without_pm<Iter1, Iter2>(
    s1: Iter1,
    mut len1: usize,
    s2: Iter2,
    mut len2: usize,
    score_cutoff: f64,
) -> Option<f64>
where
    Iter1: Iterator + Clone,
    Iter2: Iterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    let len1_orig = len1;
    let len2_orig = len2;

    if score_cutoff > 1.0 {
        return None;
    }

    if len1_orig == 0 && len2_orig == 0 {
        return Some(1.0);
    }

    // filter out based on the length difference between the two strings
    if !length_filter(len1_orig, len2_orig, score_cutoff) {
        return if score_cutoff > 0.0 { None } else { Some(0.0) };
    }

    if len1_orig == 1 && len2_orig == 1 {
        return Some(s1.eq(s2).into());
    }

    // since jaro uses a sliding window some parts of T/P might never be in
    // range an can be removed ahead of time
    let bound;
    if len2 > len1 {
        bound = len2 / 2 - 1;

        if len2 > len1 + bound {
            len2 = len1 + bound;
        }
    } else {
        bound = len1 / 2 - 1;

        if len1 > len2 + bound {
            len1 = len2 + bound;
        }
    };
    let s1_iter_win = s1.take(len1);
    let s2_iter_win = s2.take(len2);

    // common prefix never includes Transpositions
    let mut common_chars = find_common_prefix(s1_iter_win.clone(), s2_iter_win.clone());
    let s1_iter = s1_iter_win.skip(common_chars);
    let s2_iter = s2_iter_win.skip(common_chars);
    len1 -= common_chars;
    len2 -= common_chars;
    let mut transpositions = 0_usize;

    if len1 == 0 || len2 == 0 {
        // already has correct number of common chars and transpositions
    } else if len1 <= 64 && len2 <= 64 {
        let mut pm = PatternMatchVector::default();
        pm.insert(s1_iter);

        let flagged = flag_similar_characters_word(&pm, len1, s2_iter.clone(), len2, bound);

        common_chars += flagged.count_common_chars();

        if !common_char_filter(len1_orig, len2_orig, common_chars, score_cutoff) {
            return if score_cutoff > 0.0 { None } else { Some(0.0) };
        }

        transpositions = count_transpositions_word(&pm, s2_iter, len2, &flagged);
    } else {
        let mut pm = BlockPatternMatchVector::new(len1);
        pm.insert(s1_iter);

        let flagged = flag_similar_characters_block(&pm, len1, s2_iter.clone(), len2, bound);

        let flagged_chars = flagged.count_common_chars();
        common_chars += flagged_chars;

        if !common_char_filter(len1_orig, len2_orig, common_chars, score_cutoff) {
            return if score_cutoff > 0.0 { None } else { Some(0.0) };
        }

        transpositions = count_transpositions_block(&pm, s2_iter, len2, &flagged, flagged_chars);
    }

    let sim = calculate_similarity(len1_orig, len2_orig, common_chars, transpositions);
    if sim >= score_cutoff {
        Some(sim)
    } else {
        None
    }
}

pub(crate) fn similarity_with_pm<Iter1, Iter2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    mut len1: usize,
    s2: Iter2,
    mut len2: usize,
    score_cutoff: f64,
) -> Option<f64>
where
    Iter1: Iterator + Clone,
    Iter2: Iterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    let len1_orig = len1;
    let len2_orig = len2;

    if score_cutoff > 1.0 {
        return None;
    }

    if len1_orig == 0 && len2_orig == 0 {
        return Some(1.0);
    }

    // filter out based on the length difference between the two strings
    if !length_filter(len1_orig, len2_orig, score_cutoff) {
        return if score_cutoff > 0.0 { None } else { Some(0.0) };
    }

    if len1_orig == 1 && len2_orig == 1 {
        return Some(s1.eq(s2).into());
    }

    // since jaro uses a sliding window some parts of T/P might never be in
    // range an can be removed ahead of time
    let bound;
    if len2 > len1 {
        bound = len2 / 2 - 1;

        if len2 > len1 + bound {
            len2 = len1 + bound;
        }
    } else {
        bound = len1 / 2 - 1;

        if len1 > len2 + bound {
            len1 = len2 + bound;
        }
    };
    let s2_iter = s2.take(len2);

    // common prefix never includes Transpositions
    let mut common_chars = 0_usize;
    let mut transpositions = 0_usize;

    if len1 == 0 || len2 == 0 {
        // already has correct number of common chars and transpositions
    } else if len1 <= 64 && len2 <= 64 {
        let flagged = flag_similar_characters_word(pm, len1, s2_iter.clone(), len2, bound);

        common_chars += flagged.count_common_chars();

        if !common_char_filter(len1_orig, len2_orig, common_chars, score_cutoff) {
            return if score_cutoff > 0.0 { None } else { Some(0.0) };
        }

        transpositions = count_transpositions_word(pm, s2_iter, len2, &flagged);
    } else {
        let flagged = flag_similar_characters_block(pm, len1, s2_iter.clone(), len2, bound);

        let flagged_chars = flagged.count_common_chars();
        common_chars += flagged_chars;

        if !common_char_filter(len1_orig, len2_orig, common_chars, score_cutoff) {
            return if score_cutoff > 0.0 { None } else { Some(0.0) };
        }

        transpositions = count_transpositions_block(pm, s2_iter, len2, &flagged, flagged_chars);
    }

    let sim = calculate_similarity(len1_orig, len2_orig, common_chars, transpositions);
    if sim >= score_cutoff {
        Some(sim)
    } else {
        None
    }
}

pub(crate) struct IndividualComparator;

impl Metricf64 for IndividualComparator {
    fn maximum(&self, _len1: usize, _len2: usize) -> f64 {
        1.0
    }

    fn _similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        _score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        similarity_without_pm(s1, len1, s2, len2, score_cutoff.unwrap_or(0.0))
    }
}

/// Jaro distance in the range [0.0, 1.0].
///
/// This is calculated as `1.0 - `[`similarity`].
///
pub fn distance<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
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

/// Jaro similarity in the range [1.0, 0.0].
pub fn similarity<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    ScoreCutoff: Into<Option<f64>>,
    ScoreHint: Into<Option<f64>>,
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

/// Normalized Jaro distance in the range [0.0, 1.0].
///
/// This behaves the same as `distance`, since the Jaro similarity is always
/// normalized
///
pub fn normalized_distance<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
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

/// Normalized Jaro similarity in the range [1.0, 0.0].
///
/// This behaves the same as `similarity`, since the Jaro similarity is always
/// normalized
///
pub fn normalized_similarity<Iter1, Iter2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    score_cutoff: ScoreCutoff,
    score_hint: ScoreHint,
) -> Option<f64>
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
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

/// `One x Many` comparisons using the Jaro similarity
pub struct BatchComparator<Elem1> {
    s1: Vec<Elem1>,
    pm: BlockPatternMatchVector,
}

impl<CharT> Metricf64 for BatchComparator<CharT> {
    fn maximum(&self, _len1: usize, _len2: usize) -> f64 {
        1.0
    }

    fn _similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        _score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        similarity_with_pm(&self.pm, s1, len1, s2, len2, score_cutoff.unwrap_or(0.0))
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
    pub fn normalized_distance<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
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
    pub fn normalized_similarity<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
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
    pub fn distance<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
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
    pub fn similarity<Iter2, ScoreCutoff, ScoreHint>(
        &self,
        s2: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>>,
        ScoreHint: Into<Option<f64>>,
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

    fn _test_distance<Iter1, Iter2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
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

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn _test_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = _test_distance(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = _test_distance(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    fn _test_similarity<Iter1, Iter2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
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

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn _test_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = _test_similarity(
            s1.chars(),
            s2.chars(),
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = _test_similarity(s1.bytes(), s2.bytes(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    #[test]
    fn test_no_cutoff() {
        assert_delta!(
            Some(0.455556),
            _test_similarity_ascii("james", "robert", None, None),
            0.0001
        );
        assert_delta!(
            Some(1.0 - 0.455556),
            _test_distance_ascii("james", "robert", None, None),
            0.0001
        );
    }

    #[test]
    fn test_flag_chars() {
        let names = [
            "james",
            "robert",
            "john",
            "michael",
            "william",
            "david",
            "joseph",
            "thomas",
            "charles",
            "mary",
            "patricia",
            "jennifer",
            "linda",
            "elizabeth",
            "barbara",
            "susan",
            "jessica",
            "sarah",
            "karen",
            "",
        ];

        let score_cutoffs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1];

        let scores = [
            1.0, 0.455556, 0.483333, 0.561905, 0.0, 0.466667, 0.588889, 0.577778, 0.67619,
            0.483333, 0.441667, 0.55, 0.0, 0.374074, 0.447619, 0.0, 0.67619, 0.466667, 0.6, 0.0,
            0.455556, 1.0, 0.472222, 0.436508, 0.0, 0.0, 0.555556, 0.444444, 0.373016, 0.472222,
            0.361111, 0.527778, 0.0, 0.5, 0.531746, 0.0, 0.436508, 0.455556, 0.577778, 0.0,
            0.483333, 0.472222, 1.0, 0.464286, 0.0, 0.0, 0.611111, 0.444444, 0.464286, 0.0, 0.0,
            0.583333, 0.483333, 0.0, 0.0, 0.483333, 0.464286, 0.0, 0.483333, 0.0, 0.561905,
            0.436508, 0.464286, 1.0, 0.52381, 0.447619, 0.373016, 0.539683, 0.742857, 0.464286,
            0.490079, 0.511905, 0.561905, 0.587302, 0.428571, 0.447619, 0.428571, 0.395238,
            0.447619, 0.0, 0.0, 0.0, 0.0, 0.52381, 1.0, 0.447619, 0.0, 0.436508, 0.428571, 0.0,
            0.60119, 0.422619, 0.565079, 0.47619, 0.428571, 0.447619, 0.52381, 0.447619, 0.0, 0.0,
            0.466667, 0.0, 0.0, 0.447619, 0.447619, 1.0, 0.0, 0.0, 0.447619, 0.483333, 0.55,
            0.441667, 0.466667, 0.374074, 0.447619, 0.0, 0.447619, 0.466667, 0.466667, 0.0,
            0.588889, 0.555556, 0.611111, 0.373016, 0.0, 0.0, 1.0, 0.444444, 0.436508, 0.0, 0.0,
            0.527778, 0.0, 0.518519, 0.0, 0.455556, 0.531746, 0.577778, 0.455556, 0.0, 0.577778,
            0.444444, 0.444444, 0.539683, 0.436508, 0.0, 0.444444, 1.0, 0.642857, 0.0, 0.361111,
            0.0, 0.455556, 0.425926, 0.436508, 0.455556, 0.373016, 0.455556, 0.0, 0.0, 0.67619,
            0.373016, 0.464286, 0.742857, 0.428571, 0.447619, 0.436508, 0.642857, 1.0, 0.595238,
            0.511905, 0.422619, 0.447619, 0.47619, 0.52381, 0.447619, 0.0, 0.561905, 0.67619, 0.0,
            0.483333, 0.472222, 0.0, 0.464286, 0.0, 0.483333, 0.0, 0.0, 0.595238, 1.0, 0.583333,
            0.0, 0.0, 0.453704, 0.595238, 0.0, 0.0, 0.633333, 0.633333, 0.0, 0.441667, 0.361111,
            0.0, 0.490079, 0.60119, 0.55, 0.0, 0.361111, 0.511905, 0.583333, 1.0, 0.416667,
            0.383333, 0.324074, 0.60119, 0.441667, 0.60119, 0.55, 0.55, 0.0, 0.55, 0.527778,
            0.583333, 0.511905, 0.422619, 0.441667, 0.527778, 0.0, 0.422619, 0.0, 0.416667, 1.0,
            0.383333, 0.569444, 0.422619, 0.441667, 0.60119, 0.0, 0.55, 0.0, 0.0, 0.0, 0.483333,
            0.561905, 0.565079, 0.466667, 0.0, 0.455556, 0.447619, 0.0, 0.383333, 0.383333, 1.0,
            0.644444, 0.447619, 0.466667, 0.447619, 0.466667, 0.0, 0.0, 0.374074, 0.5, 0.0,
            0.587302, 0.47619, 0.374074, 0.518519, 0.425926, 0.47619, 0.453704, 0.324074, 0.569444,
            0.644444, 1.0, 0.502646, 0.437037, 0.587302, 0.437037, 0.374074, 0.0, 0.447619,
            0.531746, 0.0, 0.428571, 0.428571, 0.447619, 0.0, 0.436508, 0.52381, 0.595238, 0.60119,
            0.422619, 0.447619, 0.502646, 1.0, 0.447619, 0.428571, 0.67619, 0.561905, 0.0, 0.0,
            0.0, 0.483333, 0.447619, 0.447619, 0.0, 0.455556, 0.455556, 0.447619, 0.0, 0.441667,
            0.441667, 0.466667, 0.437037, 0.447619, 1.0, 0.561905, 0.6, 0.466667, 0.0, 0.67619,
            0.436508, 0.464286, 0.428571, 0.52381, 0.447619, 0.531746, 0.373016, 0.0, 0.0, 0.60119,
            0.60119, 0.447619, 0.587302, 0.428571, 0.561905, 1.0, 0.447619, 0.447619, 0.0,
            0.466667, 0.455556, 0.0, 0.395238, 0.447619, 0.466667, 0.577778, 0.455556, 0.561905,
            0.633333, 0.55, 0.0, 0.466667, 0.437037, 0.67619, 0.6, 0.447619, 1.0, 0.6, 0.0, 0.6,
            0.577778, 0.483333, 0.447619, 0.0, 0.466667, 0.455556, 0.0, 0.67619, 0.633333, 0.55,
            0.55, 0.0, 0.374074, 0.561905, 0.466667, 0.447619, 0.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ];

        for score_cutoff in score_cutoffs {
            for (i, name1) in names.iter().enumerate() {
                for (j, name2) in names.iter().enumerate() {
                    let score = scores[i * names.len() + j];
                    let expected_sim = if score_cutoff <= score {
                        Some(score)
                    } else {
                        None
                    };
                    let expected_dist = expected_sim.map(|s| 1.0 - s);

                    let sim = _test_similarity_ascii(name1, name2, score_cutoff, None);
                    let dist = _test_distance_ascii(name1, name2, 1.0 - score_cutoff, None);
                    assert_delta!(expected_sim, sim, 0.0001);
                    assert_delta!(expected_dist, dist, 0.0001);
                }
            }
        }
    }
}
