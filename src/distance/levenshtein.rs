use crate::details::common::{remove_common_affix, HashableChar};
use crate::details::distance::MetricUsize;
use crate::details::growing_hashmap::{GrowingHashmap, HybridGrowingHashmap};
use crate::details::intrinsics::{ceil_div_usize, shr64};
use crate::details::matrix::ShiftedBitMatrix;
use crate::details::pattern_match_vector::{
    BitVectorInterface, BitvectorHashmap, BlockPatternMatchVector, PatternMatchVector,
};
use crate::distance::indel;
use std::cmp::{max, min};
use std::mem;

#[derive(Clone, Copy)]
pub struct WeightTable {
    pub insert_cost: usize,
    pub delete_cost: usize,
    pub replace_cost: usize,
}

impl Default for WeightTable {
    fn default() -> Self {
        Self {
            insert_cost: 1,
            delete_cost: 1,
            replace_cost: 1,
        }
    }
}

#[derive(Clone)]
struct LevenshteinRow {
    vp: u64,
    vn: u64,
}

impl Default for LevenshteinRow {
    fn default() -> Self {
        Self { vp: !0_u64, vn: 0 }
    }
}

#[derive(Default)]
struct ResultMatrix {
    vp: ShiftedBitMatrix<u64>,
    vn: ShiftedBitMatrix<u64>,
}

#[derive(Default)]
struct ResultRow {
    first_block: usize,
    last_block: usize,
    prev_score: usize,
    vecs: Vec<LevenshteinRow>,
}

struct DistanceResult<const RECORD_MATRIX: usize, const RECORD_BIT_ROW: usize> {
    record_matrix: [ResultMatrix; RECORD_MATRIX],
    bit_row: [ResultRow; RECORD_BIT_ROW],
    dist: Option<usize>,
}

impl Default for DistanceResult<0, 0> {
    fn default() -> Self {
        Self {
            record_matrix: [],
            bit_row: [],
            dist: None,
        }
    }
}

impl Default for DistanceResult<1, 0> {
    fn default() -> Self {
        Self {
            record_matrix: [ResultMatrix::default()],
            bit_row: [],
            dist: None,
        }
    }
}

impl Default for DistanceResult<0, 1> {
    fn default() -> Self {
        Self {
            record_matrix: [],
            bit_row: [ResultRow::default()],
            dist: None,
        }
    }
}

fn generalized_wagner_fischer<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    _len2: usize,
    weights: &WeightTable,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2>,
    Elem2: PartialEq<Elem1>,
{
    let cache_size = len1 + 1;
    let mut cache: Vec<usize> = (0..cache_size).map(|x| x * weights.delete_cost).collect();

    for ch2 in s2 {
        let mut cache_iter = cache.iter_mut().peekable();
        let mut cur_cache: &mut usize = cache_iter
            .next()
            .expect("cache always has at least one element");
        let mut temp = *cur_cache;
        *cur_cache += weights.insert_cost;

        for ch1 in s1.clone() {
            if ch1 != ch2 {
                temp = min(
                    *cur_cache + weights.delete_cost,
                    temp + weights.replace_cost,
                );
                temp = min(
                    temp,
                    **cache_iter
                        .peek()
                        .expect("cache has len1 + 1 elements, so this should always exist")
                        + weights.insert_cost,
                );
            }

            cur_cache = cache_iter
                .next()
                .expect("cache has len1 + 1 elements, so this should always exist");
            mem::swap(cur_cache, &mut temp);
        }
    }

    let dist = *cache.last().expect("cache always has at least one element");
    if dist <= score_cutoff {
        Some(dist)
    } else {
        None
    }
}

/**
 * @brief calculates the maximum possible Levenshtein distance based on
 * string lengths and weights
 */
fn _maximum(len1: usize, len2: usize, weights: &WeightTable) -> usize {
    let max_dist = len1 * weights.delete_cost + len2 * weights.insert_cost;

    if len1 >= len2 {
        min(
            max_dist,
            len2 * weights.replace_cost + (len1 - len2) * weights.delete_cost,
        )
    } else {
        min(
            max_dist,
            len1 * weights.replace_cost + (len2 - len1) * weights.insert_cost,
        )
    }
}

fn _min_distance(len1: usize, len2: usize, weights: &WeightTable) -> usize {
    max(
        (len1 as isize - len2 as isize) * weights.delete_cost as isize,
        (len2 as isize - len1 as isize) * weights.insert_cost as isize,
    ) as usize
}

fn generalized_distance<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    weights: &WeightTable,
    score_cutoff: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
{
    let min_edits = _min_distance(len1, len2, weights);
    if min_edits > score_cutoff {
        return None;
    }

    // common affix does not effect Levenshtein distance
    let affix = remove_common_affix(s1, len1, s2, len2);

    generalized_wagner_fischer(
        affix.s1,
        affix.len1,
        affix.s2,
        affix.len2,
        weights,
        score_cutoff,
    )
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
///   01 = DELETE, 10 = INSERT, 11 = SUBSTITUTE
///
/// For example, 3F -> 0b111111 means three substitutions
static LEVENSHTEIN_MBLEVEN2018_MATRIX: [[u8; 7]; 9] = [
    // max edit distance 1
    [0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 0
    [0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 1
    // max edit distance 2
    [0x0F, 0x09, 0x06, 0x00, 0x00, 0x00, 0x00], // len_diff 0
    [0x0D, 0x07, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 1
    [0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 2
    // max edit distance 3
    [0x3F, 0x27, 0x2D, 0x39, 0x36, 0x1E, 0x1B], // len_diff 0
    [0x3D, 0x37, 0x1F, 0x25, 0x19, 0x16, 0x00], // len_diff 1
    [0x35, 0x1D, 0x17, 0x00, 0x00, 0x00, 0x00], // len_diff 2
    [0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], // len_diff 3
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
    // todo add assert that first + last character are different
    // this is used below, but would break if common affix is not removed

    if len1 < len2 {
        return mbleven2018(s2, len2, s1, len1, score_cutoff);
    }

    let len_diff = len1 - len2;

    if score_cutoff == 1 {
        return if len_diff == 1 || len1 != 1 {
            None
        } else {
            Some(1)
        };
    }

    let ops_index = (score_cutoff + score_cutoff * score_cutoff) / 2 + len_diff - 1;
    let possible_ops = &LEVENSHTEIN_MBLEVEN2018_MATRIX[ops_index];
    let mut dist = score_cutoff + 1;

    for &ops_ in possible_ops {
        let mut ops = ops_;
        let mut iter_s1 = s1.clone();
        let mut iter_s2 = s2.clone();
        let mut cur_dist = 0;

        let mut cur1 = iter_s1.next();
        let mut cur2 = iter_s2.next();

        if ops == 0 {
            break;
        }

        loop {
            match (&cur1, &cur2) {
                (Some(ch1), Some(ch2)) => {
                    if ch1 == ch2 {
                        cur1 = iter_s1.next();
                        cur2 = iter_s2.next();
                    } else {
                        cur_dist += 1;
                        if ops == 0 {
                            break;
                        }
                        if (ops & 1) != 0 {
                            cur1 = iter_s1.next();
                        }
                        if (ops & 2) != 0 {
                            cur2 = iter_s2.next();
                        }

                        ops >>= 2;
                    }
                }
                (Some(_), None) => {
                    cur_dist += 1;
                    cur1 = iter_s1.next();
                }
                (None, Some(_)) => {
                    cur_dist += 1;
                    cur2 = iter_s2.next();
                }

                (None, None) => break,
            }
        }

        cur_dist += iter_s1.count() + iter_s2.count();
        dist = min(dist, cur_dist);
    }

    if dist <= score_cutoff {
        Some(dist)
    } else {
        None
    }
}

/// Bitparallel implementation of the Levenshtein distance.
///
/// This implementation requires the first string to have a length <= 64.
/// The algorithm used is described @cite `hyrro_2002` and has a time complexity
/// of O(N). Comments and variable names in the implementation follow the
/// paper. This implementation is used internally when the strings are short enough
fn hyrroe2003<
    const RECORD_MATRIX: usize,
    const RECORD_BIT_ROW: usize,
    PmVec,
    Iter1,
    Iter2,
    Elem1,
    Elem2,
>(
    pm: &PmVec,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> DistanceResult<RECORD_MATRIX, RECORD_BIT_ROW>
where
    Iter1: Iterator<Item = Elem1>,
    Iter2: Iterator<Item = Elem2>,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    PmVec: BitVectorInterface,
    DistanceResult<RECORD_MATRIX, RECORD_BIT_ROW>: Default,
{
    debug_assert!(len1 != 0);

    // VP is set to 1^m. Shifting by bitwidth would be undefined behavior
    let mut vp: u64 = !0_u64;
    let mut vn: u64 = 0;

    let mut dist = len1;
    let mut res = DistanceResult::<RECORD_MATRIX, RECORD_BIT_ROW>::default();
    if RECORD_MATRIX == 1 {
        res.record_matrix[0].vp = ShiftedBitMatrix::<u64>::new(len2, 1, !0_u64);
        res.record_matrix[0].vn = ShiftedBitMatrix::<u64>::new(len2, 1, 0);
    }
    // mask used when computing D[m,j] in the paper 10^(m-1)
    let mask: u64 = 1_u64 << (len1 - 1);

    for (i, ch2) in s2.enumerate() {
        let pm_j = pm.get(0, ch2);
        let x = pm_j;
        let d0 = ((x & vp).wrapping_add(vp) ^ vp) | x | vn;

        // Step 2: Computing HP and HN
        let mut hp = vn | !(d0 | vp);
        let mut hn = d0 & vp;

        // Step 3: Computing the value D[m,j]
        dist += usize::from(hp & mask != 0);
        dist -= usize::from(hn & mask != 0);

        // Step 4: Computing Vp and VN
        hp = (hp << 1) | 1;
        hn <<= 1;

        vp = hn | !(d0 | hp);
        vn = hp & d0;

        if RECORD_MATRIX == 1 {
            *res.record_matrix[0].vp.get_mut(i, 0) = vp;
            *res.record_matrix[0].vn.get_mut(i, 0) = vn;
        }
    }

    res.dist = if dist <= score_cutoff {
        Some(dist)
    } else {
        None
    };

    if RECORD_BIT_ROW == 1 {
        res.bit_row[0].first_block = 0;
        res.bit_row[0].last_block = 0;

        res.bit_row[0].prev_score = len2;
        res.bit_row[0].vecs.push(LevenshteinRow { vp, vn });
    }

    res
}

fn hyrroe2003_small_band_with_pm<PmVec, Iter1, Iter2, Elem1, Elem2>(
    pm: &PmVec,
    _s1: Iter1,
    len1: usize,
    mut s2: Iter2,
    len2: usize,
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
    let mut vp: u64 = !0_u64 << (64 - score_cutoff - 1);
    let mut vn: u64 = 0;

    let words = pm.size();
    let mut curr_dist = score_cutoff;
    let diagonal_mask = 1_u64 << 63;
    let mut horizontal_mask = 1_u64 << 62;
    let mut start_pos = (score_cutoff as isize) + 1 - 64;

    // score can decrease along the horizontal, but not along the diagonal
    let break_score =
        (score_cutoff as isize + len2 as isize - (len1 as isize - score_cutoff as isize)) as usize;

    if len1 > score_cutoff {
        for ch2 in s2.by_ref().take(len1 - score_cutoff) {
            // Step 1: Computing D0
            let mut pm_j: u64;
            if start_pos < 0 {
                pm_j = pm.get(0, ch2) << (-start_pos);
            } else {
                let word = start_pos as usize / 64;
                let word_pos = start_pos as usize % 64;

                pm_j = pm.get(word, ch2) >> word_pos;

                if word + 1 < words && word_pos != 0 {
                    pm_j |= pm.get(word + 1, ch2) << (64 - word_pos);
                }
            }

            let x = pm_j;
            let d0 = ((x & vp).wrapping_add(vp) ^ vp) | x | vn;

            // Step 2: Computing HP and HN
            let hp = vn | !(d0 | vp);
            let hn = d0 & vp;

            // Step 3: Computing the value D[m,j]
            curr_dist += usize::from(d0 & diagonal_mask == 0);

            if curr_dist > break_score {
                return None;
            }

            // Step 4: Computing Vp and VN
            vp = hn | !((d0 >> 1) | hp);
            vn = (d0 >> 1) & hp;

            start_pos += 1;
        }
    }

    for ch2 in s2 {
        // Step 1: Computing D0
        let mut pm_j: u64;
        if start_pos < 0 {
            pm_j = pm.get(0, ch2) << (-start_pos);
        } else {
            let word = start_pos as usize / 64;
            let word_pos = start_pos as usize % 64;

            pm_j = pm.get(word, ch2) >> word_pos;

            if word + 1 < words && word_pos != 0 {
                pm_j |= pm.get(word + 1, ch2) << (64 - word_pos);
            }
        }

        let x = pm_j;
        let d0 = ((x & vp).wrapping_add(vp) ^ vp) | x | vn;

        // Step 2: Computing HP and HN
        let hp = vn | !(d0 | vp);
        let hn = d0 & vp;

        // Step 3: Computing the value D[m,j]
        curr_dist += usize::from(hp & horizontal_mask != 0);
        curr_dist -= usize::from(hn & horizontal_mask != 0);
        horizontal_mask >>= 1;

        if curr_dist > break_score {
            return None;
        }

        // Step 4: Computing Vp and VN
        vp = hn | !((d0 >> 1) | hp);
        vn = (d0 >> 1) & hp;

        start_pos += 1;
    }

    if curr_dist <= score_cutoff {
        Some(curr_dist)
    } else {
        None
    }
}

fn hyrroe2003_small_band_without_pm<const RECORD_MATRIX: usize, Iter1, Iter2, Elem1, Elem2>(
    mut s1: Iter1,
    len1: usize,
    mut s2: Iter2,
    len2: usize,
    score_cutoff: usize,
) -> DistanceResult<RECORD_MATRIX, 0>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar,
    Elem2: PartialEq<Elem1> + HashableChar,
    DistanceResult<RECORD_MATRIX, 0>: Default,
{
    debug_assert!(score_cutoff <= len1);
    debug_assert!(score_cutoff <= len2);
    debug_assert!(len2 >= len1 - score_cutoff);

    // VP is set to 1^m. Shifting by bitwidth would be undefined behavior
    let mut vp: u64 = !0_u64 << (64 - score_cutoff - 1);
    let mut vn: u64 = 0;

    let mut dist = score_cutoff;
    let mut res = DistanceResult::<RECORD_MATRIX, 0>::default();
    if RECORD_MATRIX == 1 {
        res.record_matrix[0].vp = ShiftedBitMatrix::<u64>::new(len2, 1, !0_u64);
        res.record_matrix[0].vn = ShiftedBitMatrix::<u64>::new(len2, 1, 0);

        let start_offset = score_cutoff as isize + 2 - 64;
        for i in 0..len2 {
            res.record_matrix[0]
                .vp
                .set_offset(i, start_offset + i as isize);
            res.record_matrix[0]
                .vn
                .set_offset(i, start_offset + i as isize);
        }
    }

    let diagonal_mask = 1_u64 << 63;
    let mut horizontal_mask = 1_u64 << 62;

    // score can decrease along the horizontal, but not along the diagonal
    let break_score =
        (score_cutoff as isize + len2 as isize - (len1 as isize - score_cutoff as isize)) as usize;
    // rust fails to remove the array memcpy on return
    // let pm = HybridGrowingHashmap::<(isize, u64)>::new();
    let mut pm = HybridGrowingHashmap::<(isize, u64)> {
        map_unsigned: GrowingHashmap::default(),
        map_signed: GrowingHashmap::default(),
        extended_ascii: [Default::default(); 256],
    };

    let mut i = 0 - score_cutoff as isize;
    for ch1 in s1.by_ref().take(score_cutoff) {
        let item = pm.get_mut(ch1);
        item.1 = shr64(item.1, (i - item.0) as usize) | (1_u64 << 63);
        item.0 = i;
        i += 1;
    }

    // Searching
    for (ch1, ch2) in (&mut s1).zip(&mut s2).take(len1 - score_cutoff) {
        // Step 1: Computing D0
        // update bitmasks online
        {
            let item = pm.get_mut(ch1);
            item.1 = shr64(item.1, (i - item.0) as usize) | (1_u64 << 63);
            item.0 = i;
        }
        let pm_j: u64;
        {
            let item = pm.get(ch2);
            pm_j = shr64(item.1, (i - item.0) as usize);
        }

        let x = pm_j;
        let d0 = ((x & vp).wrapping_add(vp) ^ vp) | x | vn;

        // Step 2: Computing HP and HN
        let hp = vn | !(d0 | vp);
        let hn = d0 & vp;

        // Step 3: Computing the value D[m,j]
        dist += usize::from(d0 & diagonal_mask == 0);

        if dist > break_score {
            res.dist = None;
            return res;
        }

        // Step 4: Computing Vp and VN
        vp = hn | !((d0 >> 1) | hp);
        vn = (d0 >> 1) & hp;

        if RECORD_MATRIX == 1 {
            *res.record_matrix[0].vp.get_mut(i as usize, 0) = vp;
            *res.record_matrix[0].vn.get_mut(i as usize, 0) = vn;
        }
        i += 1;
    }

    for ch2 in s2 {
        if let Some(ch1) = s1.next() {
            let item = pm.get_mut(ch1);
            item.1 = shr64(item.1, (i - item.0) as usize) | (1_u64 << 63);
            item.0 = i;
        }
        let pm_j: u64;
        {
            let item = pm.get(ch2);
            pm_j = shr64(item.1, (i - item.0) as usize);
        }

        let x = pm_j;
        let d0 = ((x & vp).wrapping_add(vp) ^ vp) | x | vn;

        // Step 2: Computing HP and HN
        let hp = vn | !(d0 | vp);
        let hn = d0 & vp;

        // Step 3: Computing the value D[m,j]
        dist += usize::from(hp & horizontal_mask != 0);
        dist -= usize::from(hn & horizontal_mask != 0);
        horizontal_mask >>= 1;

        if dist > break_score {
            res.dist = None;
            return res;
        }

        // Step 4: Computing Vp and VN
        vp = hn | !((d0 >> 1) | hp);
        vn = (d0 >> 1) & hp;

        if RECORD_MATRIX == 1 {
            *res.record_matrix[0].vp.get_mut(i as usize, 0) = vp;
            *res.record_matrix[0].vn.get_mut(i as usize, 0) = vn;
        }
        i += 1;
    }

    res.dist = if dist <= score_cutoff {
        Some(dist)
    } else {
        None
    };
    res
}

#[allow(clippy::too_many_lines)]
fn hyrroe2003_block<
    const RECORD_MATRIX: usize,
    const RECORD_BIT_ROW: usize,
    Iter1,
    Iter2,
    Elem1,
    Elem2,
>(
    pm: &BlockPatternMatchVector,
    _s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    mut score_cutoff: usize,
    stop_row: isize,
) -> DistanceResult<RECORD_MATRIX, RECORD_BIT_ROW>
where
    Iter1: Iterator<Item = Elem1> + Clone,
    Iter2: Iterator<Item = Elem2> + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
    DistanceResult<RECORD_MATRIX, RECORD_BIT_ROW>: Default,
{
    let mut res: DistanceResult<RECORD_MATRIX, RECORD_BIT_ROW> = DistanceResult::default();
    if score_cutoff < len1.abs_diff(len2) {
        res.dist = None;
        return res;
    }

    let word_size = 64;
    let words = pm.size();
    let mut vecs = vec![LevenshteinRow::default(); words];
    let mut scores: Vec<usize> = (0..words).map(|x| (x + 1) * word_size).collect();
    scores[words - 1] = len1;
    let last = 1_u64 << ((len1 - 1) % word_size);

    if RECORD_MATRIX == 1 {
        let full_band = min(len1, 2 * score_cutoff + 1);
        let full_band_words = min(words, full_band / word_size + 2);
        res.record_matrix[0].vp = ShiftedBitMatrix::<u64>::new(len2, full_band_words, !0_u64);
        res.record_matrix[0].vn = ShiftedBitMatrix::<u64>::new(len2, full_band_words, 0);
    }
    if RECORD_BIT_ROW == 1 {
        res.bit_row[0].first_block = 0;
        res.bit_row[0].last_block = 0;
        res.bit_row[0].prev_score = 0;
    }

    score_cutoff = min(score_cutoff, max(len1, len2));
    // first_block is the index of the first block in Ukkonen band.
    let mut first_block: usize = 0;
    // last_block is the index of the last block in Ukkonen band.
    let mut last_block = min(
        words,
        ceil_div_usize(
            min(score_cutoff, (score_cutoff + len1 - len2) / 2) + 1,
            word_size,
        ),
    ) - 1;

    // Searching
    for (row, ch2) in s2.enumerate() {
        let mut hp_carry: bool = true;
        let mut hn_carry: bool = false;

        if RECORD_MATRIX == 1 {
            res.record_matrix[0]
                .vp
                .set_offset(row, (first_block * word_size) as isize);
            res.record_matrix[0]
                .vn
                .set_offset(row, (first_block * word_size) as isize);
        }

        // rust can't mutably borrow + borrow at the same time and the closure would borrow on creation.
        // so pass everything we need to borrow in explicitly ...
        let mut advance_block = |word: usize,
                                 vecs_: &mut Vec<LevenshteinRow>,
                                 hp_carry_: &mut bool,
                                 hn_carry_: &mut bool| {
            // Step 1: Computing D0
            let pm_j = pm.get(word, ch2);
            let vn = vecs_[word].vn;
            let vp = vecs_[word].vp;

            let x = pm_j | u64::from(*hn_carry_);
            let d0 = ((x & vp).wrapping_add(vp) ^ vp) | x | vn;

            // Step 2: Computing HP and HN
            let mut hp = vn | !(d0 | vp);
            let mut hn = d0 & vp;

            let hp_carry_temp = *hp_carry_;
            let hn_carry_temp = *hn_carry_;
            if word < words - 1 {
                *hp_carry_ = (hp >> 63) != 0;
                *hn_carry_ = (hn >> 63) != 0;
            } else {
                *hp_carry_ = (hp & last) != 0;
                *hn_carry_ = (hn & last) != 0;
            }

            // Step 4: Computing Vp and VN
            hp = (hp << 1) | u64::from(hp_carry_temp);
            hn = (hn << 1) | u64::from(hn_carry_temp);

            vecs_[word].vp = hn | !(d0 | hp);
            vecs_[word].vn = hp & d0;

            if RECORD_MATRIX == 1 {
                *res.record_matrix[0].vp.get_mut(row, word - first_block) = vecs_[word].vp;
                *res.record_matrix[0].vn.get_mut(row, word - first_block) = vecs_[word].vn;
            }
        };

        let get_row_num = |word: usize| -> usize {
            if word + 1 == words {
                len1 - 1
            } else {
                (word + 1) * word_size - 1
            }
        };

        for (word, score) in scores
            .iter_mut()
            .enumerate()
            .take(last_block + 1)
            .skip(first_block)
        {
            // Step 3: Computing the value D[m,j]
            advance_block(word, &mut vecs, &mut hp_carry, &mut hn_carry);
            *score += usize::from(hp_carry);
            *score -= usize::from(hn_carry);
        }

        score_cutoff = min(
            score_cutoff as isize,
            scores[last_block] as isize
                + max(
                    len2 as isize - row as isize - 1,
                    len1 as isize - ((1 + last_block) * word_size - 1) as isize - 1,
                ),
        ) as usize;

        //---------- Adjust number of blocks according to Ukkonen ----------
        // todo on the last word instead of word_size often s1.size() % 64 should be used

        // Band adjustment: last_block
        // If block is not beneath band, calculate next block. Only next because others are certainly beneath
        // band.
        if last_block + 1 < words
            && get_row_num(last_block) as isize
                <= score_cutoff as isize + 2 * word_size as isize + row as isize + len1 as isize
                    - scores[last_block] as isize
                    - 2
                    - len2 as isize
        {
            last_block += 1;
            vecs[last_block].vp = !0_u64;
            vecs[last_block].vn = 0;

            let chars_in_block = if last_block + 1 == words {
                (len1 - 1) % word_size + 1
            } else {
                64
            };
            scores[last_block] = scores[last_block - 1] + chars_in_block - usize::from(hp_carry)
                + usize::from(hn_carry);

            advance_block(last_block, &mut vecs, &mut hp_carry, &mut hn_carry);
            scores[last_block] += usize::from(hp_carry);
            scores[last_block] -= usize::from(hn_carry);
        }

        while last_block >= first_block {
            // in band if score <= k where score >= score_last - word_size + 1
            let in_band_cond1 = scores[last_block] < score_cutoff + word_size;

            // in band if row <= max - score - len2 + len1 + i
            // if the condition is met for the first cell in the block, it
            // is met for all other cells in the blocks as well
            //
            // this uses a more loose condition similar to edlib:
            // https://github.com/Martinsos/edlib
            let in_band_cond2 = get_row_num(last_block) as isize
                <= score_cutoff as isize
                    + 2 * word_size as isize
                    + row as isize
                    + len1 as isize
                    + 1
                    - scores[last_block] as isize
                    - 2
                    - len2 as isize;

            if in_band_cond1 && in_band_cond2 {
                break;
            }
            last_block -= 1;
        }

        // Band adjustment: first_block
        while first_block <= last_block {
            // in band if score <= k where score >= score_last - word_size + 1
            let in_band_cond1 = scores[first_block] < score_cutoff + word_size;

            // in band if row >= score - max - len2 + len1 + i
            // if this condition is met for the last cell in the block, it
            // is met for all other cells in the blocks as well
            let in_band_cond2 = get_row_num(first_block) as isize
                >= scores[first_block] as isize + len1 as isize + row as isize
                    - score_cutoff as isize
                    - len2 as isize;

            if in_band_cond1 && in_band_cond2 {
                break;
            }
            first_block += 1;
        }

        // distance is larger than max, so band stops to exist
        if last_block < first_block {
            res.dist = None;
            return res;
        }

        if RECORD_BIT_ROW == 1 && row as isize == stop_row {
            if first_block == 0 {
                res.bit_row[0].prev_score = row + 1;
            } else {
                // count backwards to find score at last position in previous block
                let relevant_bits = min((first_block + 1) * 64, len1) % 64;
                let mut mask = !0_u64;
                if relevant_bits != 0 {
                    mask >>= 64 - relevant_bits;
                }

                res.bit_row[0].prev_score = scores[first_block]
                    + (vecs[first_block].vn & mask).count_ones() as usize
                    - (vecs[first_block].vp & mask).count_ones() as usize;
            }

            res.bit_row[0].first_block = first_block;
            res.bit_row[0].last_block = last_block;
            mem::swap(&mut res.bit_row[0].vecs, &mut vecs);
            // unknown so set to None
            res.dist = None;
            return res;
        }
    }

    let dist = scores[words - 1];
    res.dist = if dist <= score_cutoff {
        Some(dist)
    } else {
        None
    };
    res
}

fn uniform_distance_with_pm<Iter1, Iter2, Elem1, Elem2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    mut score_cutoff: usize,
    mut score_hint: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    // upper bound
    score_cutoff = min(score_cutoff, max(len1, len2));
    score_hint = max(score_hint, 31);

    // when no differences are allowed a direct comparision is sufficient
    if score_cutoff == 0 {
        return if s1.into_iter().eq(s2) { Some(0) } else { None };
    }

    if score_cutoff < len1.abs_diff(len2) {
        return None;
    }

    // important to catch, since this causes block to be empty -> raises exception on access
    if len1 == 0 || len2 == 0 {
        return Some(len1 + len2);
    }

    // do this first, since we can not remove any affix in encoded form
    // todo actually we could at least remove the common prefix and just shift the band
    if score_cutoff >= 4 {
        // todo could safe up to 25% even without score_cutoff when ignoring irrelevant paths
        // in the upper and lower corner
        let mut full_band = min(len1, 2 * score_cutoff + 1);

        if len1 <= 64 {
            let res: DistanceResult<0, 0> = hyrroe2003(pm, s1, len1, s2, len2, score_cutoff);
            return res.dist;
        } else if full_band <= 64 {
            return hyrroe2003_small_band_with_pm(pm, s1, len1, s2, len2, score_cutoff);
        }

        while score_hint < score_cutoff {
            full_band = min(len1, 2 * score_hint + 1);

            let score = if full_band <= 64 {
                hyrroe2003_small_band_with_pm(pm, s1.clone(), len1, s2.clone(), len2, score_hint)
            } else {
                let res: DistanceResult<0, 0> =
                    hyrroe2003_block(pm, s1.clone(), len1, s2.clone(), len2, score_hint, -1);
                res.dist
            };

            if score.is_some() {
                return score;
            }

            if usize::MAX / 2 < score_hint {
                break;
            }
            score_hint *= 2;
        }

        let res: DistanceResult<0, 0> = hyrroe2003_block(pm, s1, len1, s2, len2, score_cutoff, -1);
        return res.dist;
    }

    // common affix does not effect Levenshtein distance
    let affix = remove_common_affix(s1, len1, s2, len2);

    if affix.len1 == 0 || affix.len2 == 0 {
        return Some(affix.len1 + affix.len2);
    }

    mbleven2018(affix.s1, affix.len1, affix.s2, affix.len2, score_cutoff)
}

fn uniform_distance_without_pm<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    mut score_cutoff: usize,
    mut score_hint: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    if len1 < len2 {
        return uniform_distance_without_pm(s2, len2, s1, len1, score_cutoff, score_hint);
    }

    score_cutoff = min(score_cutoff, max(len1, len2));
    score_hint = max(score_hint, 31);

    // when no differences are allowed a direct comparision is sufficient
    if score_cutoff == 0 {
        return if s1.into_iter().eq(s2) { Some(0) } else { None };
    }

    if score_cutoff < len1.abs_diff(len2) {
        return None;
    }

    // common affix does not effect Levenshtein distance
    // todo is this really the best way to remove the common affix?
    let affix = remove_common_affix(s1, len1, s2, len2);

    if affix.len1 == 0 || affix.len2 == 0 {
        return Some(affix.len1 + affix.len2);
    }

    if score_cutoff < 4 {
        return mbleven2018(affix.s1, affix.len1, affix.s2, affix.len2, score_cutoff);
    }

    // todo could safe up to 25% even without score_cutoff when ignoring irrelevant paths
    // in the upper and lower corner
    let mut full_band = min(affix.len1, 2 * score_cutoff + 1);

    /* when the short strings has less then 65 elements Hyyrös' algorithm can be used */
    if affix.len2 <= 64 {
        // rust fails to elide the copy when returning the array
        // from PatternMatchVector::new so manually inline it
        //let block = PatternMatchVector::new(s2_iter.clone());
        let mut pm = PatternMatchVector {
            map_unsigned: BitvectorHashmap::default(),
            map_signed: BitvectorHashmap::default(),
            extended_ascii: [0; 256],
        };
        pm.insert(affix.s2.clone());

        let res: DistanceResult<0, 0> = hyrroe2003(
            &pm,
            affix.s2,
            affix.len2,
            affix.s1,
            affix.len1,
            score_cutoff,
        );
        res.dist
    } else if full_band <= 64 {
        let res: DistanceResult<0, 0> = hyrroe2003_small_band_without_pm(
            affix.s1,
            affix.len1,
            affix.s2,
            affix.len2,
            score_cutoff,
        );
        res.dist
    } else {
        let mut pm = BlockPatternMatchVector::new(affix.len1);
        pm.insert(affix.s1.clone());
        while score_hint < score_cutoff {
            full_band = min(affix.len1, 2 * score_hint + 1);

            let score = if full_band <= 64 {
                hyrroe2003_small_band_with_pm(
                    &pm,
                    affix.s1.clone(),
                    affix.len1,
                    affix.s2.clone(),
                    affix.len2,
                    score_hint,
                )
            } else {
                let res: DistanceResult<0, 0> = hyrroe2003_block(
                    &pm,
                    affix.s1.clone(),
                    affix.len1,
                    affix.s2.clone(),
                    affix.len2,
                    score_hint,
                    -1,
                );
                res.dist
            };

            if score.is_some() {
                return score;
            }

            if usize::MAX / 2 < score_hint {
                break;
            }
            score_hint *= 2;
        }

        let res: DistanceResult<0, 0> = hyrroe2003_block(
            &pm,
            affix.s1.clone(),
            affix.len1,
            affix.s2.clone(),
            affix.len2,
            score_cutoff,
            -1,
        );
        res.dist
    }
}

fn _distance_without_pm<Iter1, Iter2, Elem1, Elem2>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    weights: &WeightTable,
    score_cutoff: usize,
    score_hint: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    // for very short sequences the bitparallel algorithm is not worth it
    if len1 * len2 < 90 {
        return generalized_distance(s1, len1, s2, len2, weights, score_cutoff);
    }

    if weights.insert_cost == weights.delete_cost {
        /* when insertions + deletions operations are free there can not be any edit distance */
        if weights.insert_cost == 0 {
            return Some(0);
        }

        /* uniform Levenshtein multiplied with the common factor */
        if weights.insert_cost == weights.replace_cost {
            // score_cutoff can make use of the common divisor of the three weights
            let new_score_cutoff = ceil_div_usize(score_cutoff, weights.insert_cost);
            let new_score_hint = ceil_div_usize(score_hint, weights.insert_cost);
            let mut dist =
                uniform_distance_without_pm(s1, len1, s2, len2, new_score_cutoff, new_score_hint)?;
            dist *= weights.insert_cost;
            if dist <= score_cutoff {
                return Some(dist);
            }
            return None;
        }
        /*
         * when replace_cost >= insert_cost + delete_cost no substitutions are performed
         * therefore this can be implemented as InDel distance multiplied with the common factor
         */
        else if weights.replace_cost >= weights.insert_cost + weights.delete_cost {
            // max can make use of the common divisor of the three weights
            let new_score_cutoff = ceil_div_usize(score_cutoff, weights.insert_cost);
            let new_score_hint = ceil_div_usize(score_hint, weights.insert_cost);
            let mut dist = indel::IndividualComparator {}._distance(
                s1,
                len1,
                s2,
                len2,
                Some(new_score_cutoff),
                Some(new_score_hint),
            )?;
            dist *= weights.insert_cost;
            if dist <= score_cutoff {
                return Some(dist);
            }
            return None;
        }
    }

    generalized_distance(s1, len1, s2, len2, weights, score_cutoff)
}

#[allow(clippy::too_many_arguments)]
fn _distance_with_pm<Iter1, Iter2, Elem1, Elem2>(
    pm: &BlockPatternMatchVector,
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    weights: &WeightTable,
    score_cutoff: usize,
    score_hint: usize,
) -> Option<usize>
where
    Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
    Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
    Elem1: PartialEq<Elem2> + HashableChar + Copy,
    Elem2: PartialEq<Elem1> + HashableChar + Copy,
{
    if weights.insert_cost == weights.delete_cost {
        /* when insertions + deletions operations are free there can not be any edit distance */
        if weights.insert_cost == 0 {
            return Some(0);
        }

        /* uniform Levenshtein multiplied with the common factor */
        if weights.insert_cost == weights.replace_cost {
            // score_cutoff can make use of the common divisor of the three weights
            let new_score_cutoff = ceil_div_usize(score_cutoff, weights.insert_cost);
            let new_score_hint = ceil_div_usize(score_hint, weights.insert_cost);
            let mut dist =
                uniform_distance_with_pm(pm, s1, len1, s2, len2, new_score_cutoff, new_score_hint)?;
            dist *= weights.insert_cost;
            if dist <= score_cutoff {
                return Some(dist);
            }
            return None;
        }
        /*
         * when replace_cost >= insert_cost + delete_cost no substitutions are performed
         * therefore this can be implemented as InDel distance multiplied with the common factor
         */
        else if weights.replace_cost >= weights.insert_cost + weights.delete_cost {
            // max can make use of the common divisor of the three weights
            let new_score_cutoff = ceil_div_usize(score_cutoff, weights.insert_cost);
            let mut dist = indel::distance_with_pm(pm, s1, len1, s2, len2, new_score_cutoff)?;
            dist *= weights.insert_cost;
            if dist <= score_cutoff {
                return Some(dist);
            }
            return None;
        }
    }

    generalized_distance(s1, len1, s2, len2, weights, score_cutoff)
}

struct IndividualComparator {
    weights: Option<WeightTable>,
}

impl MetricUsize for IndividualComparator {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        let weights = self.weights.unwrap_or(WeightTable {
            insert_cost: 1,
            delete_cost: 1,
            replace_cost: 1,
        });
        _maximum(len1, len2, &weights)
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let weights = self.weights.unwrap_or(WeightTable {
            insert_cost: 1,
            delete_cost: 1,
            replace_cost: 1,
        });
        _distance_without_pm(
            s1,
            len1,
            s2,
            len2,
            &weights,
            score_cutoff.unwrap_or(usize::MAX),
            score_hint.unwrap_or(usize::MAX),
        )
    }
}

pub fn distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
    s1: Iter1,
    s2: Iter2,
    weights: Option<WeightTable>,
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
    IndividualComparator { weights }._distance(
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
    weights: Option<WeightTable>,
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
    IndividualComparator { weights }._similarity(
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
    weights: Option<WeightTable>,
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
    IndividualComparator { weights }._normalized_distance(
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
    weights: Option<WeightTable>,
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
    IndividualComparator { weights }._normalized_similarity(
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
    weights: WeightTable,
}

impl<CharT> MetricUsize for BatchComparator<CharT> {
    fn maximum(&self, len1: usize, len2: usize) -> usize {
        _maximum(len1, len2, &self.weights)
    }

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        _distance_with_pm(
            &self.pm,
            s1,
            len1,
            s2,
            len2,
            &self.weights,
            score_cutoff.unwrap_or(usize::MAX),
            score_hint.unwrap_or(usize::MAX),
        )
    }
}

impl<Elem1> BatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1_: Iter1, weights_: Option<WeightTable>) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let s1_iter = s1_.into_iter();
        let s1: Vec<Elem1> = s1_iter.clone().collect();

        let weights = weights_.unwrap_or(WeightTable {
            insert_cost: 1,
            delete_cost: 1,
            replace_cost: 1,
        });

        let mut pm = BlockPatternMatchVector::new(s1.len());
        pm.insert(s1_iter);

        Self { s1, pm, weights }
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
            self.s1.iter().copied(),
            self.s1.len(),
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
            self.s1.iter().copied(),
            self.s1.len(),
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
            self.s1.iter().copied(),
            self.s1.len(),
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
    use crate::distance::example::ocr::{OCR_EXAMPLE1, OCR_EXAMPLE2};

    static EMPTY: &str = "";
    static TEST: &str = "aaaa";
    static NO_SUFFIX: &str = "aaa";
    static NO_SUFFIX2: &str = "aaab";
    static SWAPPED1: &str = "abaa";
    static SWAPPED2: &str = "baaa";
    static REPLACE_ALL: &str = "bbbb";

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

    fn _test_distance<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        weights: Option<WeightTable>,
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
            weights,
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = distance(
            s2.clone(),
            s1.clone(),
            weights,
            score_cutoff.clone(),
            score_hint.clone(),
        );

        let scorer1 = BatchComparator::new(s1.clone(), weights);
        let res3 = scorer1.distance(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone(), weights);
        let res4 = scorer2.distance(s1.clone(), score_cutoff, score_hint);

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
        assert_eq!(res1, res4);
        res1
    }

    fn _test_distance_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        weights: Option<WeightTable>,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<usize>
    where
        ScoreCutoff: Into<Option<usize>> + Clone,
        ScoreHint: Into<Option<usize>> + Clone,
    {
        let res1 = _test_distance(
            s1.chars(),
            s2.chars(),
            weights,
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = _test_distance(s1.bytes(), s2.bytes(), weights, score_cutoff, score_hint);

        assert_eq!(res1, res2);
        res1
    }

    fn _test_normalized_similarity<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
        s1_: Iter1,
        s2_: Iter2,
        weights: Option<WeightTable>,
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
            weights,
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 = normalized_similarity(
            s2.clone(),
            s1.clone(),
            weights,
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let scorer1 = BatchComparator::new(s1.clone(), weights);
        let res3 =
            scorer1.normalized_similarity(s2.clone(), score_cutoff.clone(), score_hint.clone());
        let scorer2 = BatchComparator::new(s2.clone(), weights);
        let res4 = scorer2.normalized_similarity(s1.clone(), score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        assert_delta!(res1, res3, 0.0001);
        assert_delta!(res1, res4, 0.0001);
        res1
    }

    fn _test_normalized_similarity_ascii<ScoreCutoff, ScoreHint>(
        s1: &str,
        s2: &str,
        weights: Option<WeightTable>,
        score_cutoff: ScoreCutoff,
        score_hint: ScoreHint,
    ) -> Option<f64>
    where
        ScoreCutoff: Into<Option<f64>> + Clone,
        ScoreHint: Into<Option<f64>> + Clone,
    {
        let res1 = _test_normalized_similarity(
            s1.chars(),
            s2.chars(),
            weights,
            score_cutoff.clone(),
            score_hint.clone(),
        );
        let res2 =
            _test_normalized_similarity(s1.bytes(), s2.bytes(), weights, score_cutoff, score_hint);

        assert_delta!(res1, res2, 0.0001);
        res1
    }

    /// levenshtein calculates empty sequence
    #[test]
    fn empty() {
        assert_eq!(
            Some(0),
            _test_distance_ascii(EMPTY, EMPTY, None, None, None)
        );
        assert_eq!(Some(4), _test_distance_ascii(TEST, EMPTY, None, None, None));
        assert_eq!(Some(4), _test_distance_ascii(EMPTY, TEST, None, None, None));
    }

    /// levenshtein calculates correct distances
    #[test]
    fn simple() {
        assert_eq!(Some(0), _test_distance_ascii(TEST, TEST, None, None, None));
        assert_eq!(
            Some(1),
            _test_distance_ascii(TEST, NO_SUFFIX, None, None, None)
        );
        assert_eq!(
            Some(1),
            _test_distance_ascii(TEST, NO_SUFFIX2, None, None, None)
        );
        assert_eq!(
            Some(2),
            _test_distance_ascii(SWAPPED1, SWAPPED2, None, None, None)
        );
        assert_eq!(
            Some(4),
            _test_distance_ascii(TEST, REPLACE_ALL, None, None, None)
        );

        assert_delta!(
            Some(1.0),
            _test_normalized_similarity_ascii(TEST, TEST, None, None, None),
            0.0001
        );
        assert_delta!(
            Some(0.75),
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX, None, None, None),
            0.0001
        );
        assert_delta!(
            Some(0.75),
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX2, None, None, None),
            0.0001
        );
        assert_delta!(
            Some(0.5),
            _test_normalized_similarity_ascii(SWAPPED1, SWAPPED2, None, None, None),
            0.0001
        );
        assert_delta!(
            Some(0.0),
            _test_normalized_similarity_ascii(TEST, REPLACE_ALL, None, None, None),
            0.0001
        );
    }

    /// weighted levenshtein calculates correct distances
    #[test]
    fn weighted_simple() {
        let weights = WeightTable {
            insert_cost: 1,
            delete_cost: 1,
            replace_cost: 2,
        };
        assert_eq!(
            Some(0),
            _test_distance_ascii(TEST, TEST, Some(weights), None, None)
        );
        assert_eq!(
            Some(1),
            _test_distance_ascii(TEST, NO_SUFFIX, Some(weights), None, None)
        );
        assert_eq!(
            Some(2),
            _test_distance_ascii(SWAPPED1, SWAPPED2, Some(weights), None, None)
        );
        assert_eq!(
            Some(2),
            _test_distance_ascii(TEST, NO_SUFFIX2, Some(weights), None, None)
        );
        assert_eq!(
            Some(8),
            _test_distance_ascii(TEST, REPLACE_ALL, Some(weights), None, None)
        );

        assert_delta!(
            Some(1.0),
            _test_normalized_similarity_ascii(TEST, TEST, Some(weights), None, None),
            0.0001
        );
        assert_delta!(
            Some(0.8571),
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX, Some(weights), None, None),
            0.0001
        );
        assert_delta!(
            Some(0.75),
            _test_normalized_similarity_ascii(SWAPPED1, SWAPPED2, Some(weights), None, None),
            0.0001
        );
        assert_delta!(
            Some(0.75),
            _test_normalized_similarity_ascii(TEST, NO_SUFFIX2, Some(weights), None, None),
            0.0001
        );
        assert_delta!(
            Some(0.0),
            _test_normalized_similarity_ascii(TEST, REPLACE_ALL, Some(weights), None, None),
            0.0001
        );
    }

    /// test mbleven implementation
    #[test]
    fn test_mbleven() {
        let mut a = "South Korea";
        let mut b = "North Korea";

        assert_eq!(Some(2), _test_distance_ascii(a, b, None, None, None));
        assert_eq!(Some(2), _test_distance_ascii(a, b, None, 4, None));
        assert_eq!(Some(2), _test_distance_ascii(a, b, None, 3, None));
        assert_eq!(Some(2), _test_distance_ascii(a, b, None, 2, None));
        assert_eq!(None, _test_distance_ascii(a, b, None, 1, None));
        assert_eq!(None, _test_distance_ascii(a, b, None, 0, None));

        let weights = WeightTable {
            insert_cost: 1,
            delete_cost: 1,
            replace_cost: 2,
        };
        assert_eq!(
            Some(4),
            _test_distance_ascii(a, b, Some(weights), None, None)
        );
        assert_eq!(Some(4), _test_distance_ascii(a, b, Some(weights), 4, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 3, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 2, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 1, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 0, None));

        a = "aabc";
        b = "cccd";
        assert_eq!(Some(4), _test_distance_ascii(a, b, None, None, None));
        assert_eq!(Some(4), _test_distance_ascii(a, b, None, 4, None));
        assert_eq!(None, _test_distance_ascii(a, b, None, 3, None));
        assert_eq!(None, _test_distance_ascii(a, b, None, 2, None));
        assert_eq!(None, _test_distance_ascii(a, b, None, 1, None));
        assert_eq!(None, _test_distance_ascii(a, b, None, 0, None));

        assert_eq!(
            Some(6),
            _test_distance_ascii(a, b, Some(weights), None, None)
        );
        assert_eq!(Some(6), _test_distance_ascii(a, b, Some(weights), 6, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 5, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 4, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 3, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 2, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 1, None));
        assert_eq!(None, _test_distance_ascii(a, b, Some(weights), 0, None));
    }

    /// test banded implementation
    #[test]
    fn test_banded() {
        let mut s1 = "kkkkbbbbfkkkkkkibfkkkafakkfekgkkkkkkkkkkbdbbddddddddddafkkkekkkhkk";
        let mut s2 = "khddddddddkkkkdgkdikkccccckcckkkekkkkdddddddddddafkkhckkkkkdckkkcc";
        assert_eq!(Some(36), _test_distance_ascii(s1, s2, None, None, None));
        assert_eq!(None, _test_distance_ascii(s1, s2, None, 31, None));

        s1 = "ccddcddddddddddddddddddddddddddddddddddddddddddddddddddddaaaaaaaaaaa";
        s2 = "aaaaaaaaaaaaaadddddddddbddddddddddddddddddddddddddddddddddbddddddddd";
        assert_eq!(Some(26), _test_distance_ascii(s1, s2, None, None, None));
        assert_eq!(Some(26), _test_distance_ascii(s1, s2, None, 31, None));

        s1 = "accccccccccaaaaaaaccccccccccccccccccccccccccccccacccccccccccccccccccccccccccccc\
             ccccccccccccccccccccaaaaaaaaaaaaacccccccccccccccccccccc";
        s2 = "ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc\
             ccccccccccccccccccccccccccccccccccccbcccb";
        assert_eq!(Some(24), _test_distance_ascii(s1, s2, None, None, None));
        assert_eq!(Some(24), _test_distance_ascii(s1, s2, None, 25, None));

        s1 = "miiiiiiiiiiliiiiiiibghiiaaaaaaaaaaaaaaacccfccccedddaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
              aaaaaaaaaaaaa";
        s2 = "aaaaaaajaaaaaaaabghiiaaaaaaaaaaaaaaacccfccccedddaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
              aajjdim";
        assert_eq!(Some(27), _test_distance_ascii(s1, s2, None, None, None));
        assert_eq!(Some(27), _test_distance_ascii(s1, s2, None, 27, None));

        s1 = "lllllfllllllllllllllllllllllllllllllllllllllllllllllllglllllilldcaaaaaaaaaaaaaa\
              aaaaadbbllllllllllhllllllllllllllllllllllllllgl";
        s2 = "aaaaaaaaaaaaaadbbllllllllllllllelllllllllllllllllllllllllllllllglllllilldcaaaaa\
              aaaaaaaaaaaaaadbbllllllllllllllellllllllllllllhlllllllllill";
        assert_eq!(Some(23), _test_distance_ascii(s1, s2, None, None, None));
        assert_eq!(Some(23), _test_distance_ascii(s1, s2, None, 27, None));
        assert_eq!(Some(23), _test_distance_ascii(s1, s2, None, 28, None));

        s1 = "llccacaaaaaaaaaccccccccccccccccddffaccccaccecccggggclallhcccccljif";
        s2 = "bddcbllllllbcccccccccccccccccddffccccccccebcccggggclbllhcccccljifbddcccccc";
        assert_eq!(Some(27), _test_distance_ascii(s1, s2, None, None, None));
        assert_eq!(Some(27), _test_distance_ascii(s1, s2, None, 27, None));
        assert_eq!(Some(27), _test_distance_ascii(s1, s2, None, 28, None));
    }

    #[test]
    fn test_blockwise() {
        let s1 = "a".repeat(128);
        let s2 = "b".repeat(128);
        assert_eq!(Some(128), _test_distance_ascii(&s1, &s2, None, None, None));
    }

    #[test]
    fn test_large_band() {
        assert_eq!(106514, OCR_EXAMPLE1.iter().count());
        assert_eq!(107244, OCR_EXAMPLE2.iter().count());

        assert_eq!(
            Some(5278),
            distance(OCR_EXAMPLE1.iter(), OCR_EXAMPLE2.iter(), None, None, None)
        );
        assert_eq!(
            None,
            distance(OCR_EXAMPLE1.iter(), OCR_EXAMPLE2.iter(), None, 2500, None)
        );
        assert_eq!(
            Some(5278),
            distance(OCR_EXAMPLE1.iter(), OCR_EXAMPLE2.iter(), None, None, 0)
        );
    }
}
