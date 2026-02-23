use crate::common::{
    set_decomposition, sorted_split, NoScoreCutoff, ScoreAlignment, SimilarityCutoff,
    WithScoreCutoff,
};
use crate::details::distance::MetricUsize;
use crate::details::splitted_sentence::{IsSpace, SplittedSentence};
use crate::distance::indel;
use crate::{Hash, HashableChar};
use std::collections::HashSet;

pub fn score_cutoff_to_distance(score_cutoff: f64, lensum: usize) -> usize {
    ((lensum as f64) * (1.0 - score_cutoff)).ceil() as usize
}

pub fn norm_distance(dist: usize, lensum: usize, score_cutoff: f64) -> f64 {
    let score = if lensum > 0 {
        1.0 - (dist as f64) / (lensum as f64)
    } else {
        1.0
    };

    if score >= score_cutoff {
        score
    } else {
        0.0
    }
}

#[must_use]
#[derive(Clone, Copy, Debug)]
pub struct Args<ResultType, CutoffType> {
    score_cutoff: CutoffType,
    score_hint: Option<ResultType>,
}

impl<ResultType> Default for Args<ResultType, NoScoreCutoff> {
    fn default() -> Args<ResultType, NoScoreCutoff> {
        Args {
            score_cutoff: NoScoreCutoff,
            score_hint: None,
        }
    }
}

impl<ResultType, CutoffType> Args<ResultType, CutoffType> {
    pub fn score_hint(mut self, score_hint: ResultType) -> Self {
        self.score_hint = Some(score_hint);
        self
    }

    pub fn score_cutoff(
        self,
        score_cutoff: ResultType,
    ) -> Args<ResultType, WithScoreCutoff<ResultType>> {
        Args {
            score_hint: self.score_hint,
            score_cutoff: WithScoreCutoff(score_cutoff),
        }
    }
}

/// Returns a simple ratio between two strings
pub fn ratio<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> f64
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    ratio_with_args(s1, s2, &Args::default())
}

pub fn ratio_with_args<Iter1, Iter2, CutoffType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: SimilarityCutoff<f64>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    args.score_cutoff
        .score(indel::IndividualComparator {}._normalized_similarity(
            s1_iter.clone(),
            s1_iter.count(),
            s2_iter.clone(),
            s2_iter.count(),
            args.score_cutoff.cutoff(),
            args.score_hint,
        ))
}

/// `One x Many` comparisons using `ratio`
pub struct RatioBatchComparator<Elem1> {
    scorer: indel::BatchComparator<Elem1>,
}

impl<Elem1> RatioBatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        Self {
            scorer: indel::BatchComparator::new(s1),
        }
    }

    pub fn similarity<Iter2>(&self, s2: Iter2) -> f64
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.similarity_with_args(s2, &Args::default())
    }

    pub fn similarity_with_args<Iter2, CutoffType>(
        &self,
        s2: Iter2,
        args: &Args<f64, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: SimilarityCutoff<f64>,
    {
        let s2_iter = s2.into_iter();
        args.score_cutoff
            .score(self.scorer.scorer._normalized_similarity(
                self.scorer.scorer.s1.iter().copied(),
                self.scorer.scorer.s1.len(),
                s2_iter.clone(),
                s2_iter.count(),
                args.score_cutoff.cutoff(),
                args.score_hint,
            ))
    }
}

/// Searches for the optimal alignment of the shorter string in the longer string
/// and returns the `fuzz::ratio` score for this alignment.
pub fn partial_ratio<Iter1, Iter2>(s1: Iter1, s2: Iter2, score_cutoff: f64) -> f64
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    partial_ratio_with_args(s1, s2, &Args::default().score_cutoff(score_cutoff)).unwrap_or(0.0)
}

pub fn partial_ratio_with_args<Iter1, Iter2, CutoffType>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: SimilarityCutoff<f64>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();

    let alignment = partial_ratio_alignment(
        s1_iter.clone(),
        s1_iter.count(),
        s2_iter.clone(),
        s2_iter.count(),
        args,
    );

    let score = alignment.into().map_or(0.0, |a| a.score);
    args.score_cutoff.score(score)
}

pub fn partial_ratio_alignment<Iter1, Iter2, CutoffType>(
    s1: Iter1,
    len1: usize,
    s2: Iter2,
    len2: usize,
    args: &Args<f64, CutoffType>,
) -> CutoffType::AlignmentOutput
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    CutoffType: SimilarityCutoff<f64>,
{
    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();
    let mut score_cutoff = args.score_cutoff.cutoff().unwrap_or(0.0);

    if score_cutoff > 1.0 {
        return args.score_cutoff.alignment(None);
    }

    if len1 > len2 {
        let res = partial_ratio_alignment(s2_iter.clone(), len2, s1_iter.clone(), len1, args);
        return args.score_cutoff.alignment(res.into().map(|mut alignment| {
            std::mem::swap(&mut alignment.src_start, &mut alignment.dest_start);
            std::mem::swap(&mut alignment.src_end, &mut alignment.dest_end);
            alignment
        }));
    }

    if len1 == 0 || len2 == 0 {
        return args.score_cutoff.alignment(Some(ScoreAlignment {
            score: f64::from(len1 == len2),
            src_start: 0,
            src_end: len1,
            dest_start: 0,
            dest_end: len1,
        }));
    }

    let s1_char_set = s1_iter
        .clone()
        .map(|c| c.hash_char())
        .collect::<HashSet<Hash>>();

    let indel_comp = indel::BatchComparator::new(s1_iter.clone());

    let mut res = partial_ratio_impl(
        &indel_comp,
        len1,
        s2_iter.clone(),
        len2,
        &s1_char_set,
        score_cutoff,
        args.score_hint,
    );

    if (res.score != 1.0) && (len1 == len2) {
        score_cutoff = f64::max(score_cutoff, res.score);
        let indel_comp = indel::BatchComparator::new(s2_iter.clone());
        let res2 = partial_ratio_impl(
            &indel_comp,
            len2,
            s1_iter.clone(),
            len1,
            &s1_char_set,
            score_cutoff,
            args.score_hint,
        );
        if res2.score > res.score {
            res = ScoreAlignment {
                score: res2.score,
                src_start: res2.dest_start,
                src_end: res2.dest_end,
                dest_start: res2.src_start,
                dest_end: res2.src_end,
            };
        }
    }

    let alignment = if res.score < score_cutoff {
        None
    } else {
        Some(res)
    };

    args.score_cutoff.alignment(alignment)
}

/*
implementation of partial_ratio for len(s1) <= len(s2)
*/
fn partial_ratio_impl<Elem1, Iter2>(
    comparator: &indel::BatchComparator<Elem1>,
    len1: usize,
    s2: Iter2,
    len2: usize,
    s1_char_set: &HashSet<Hash>,
    mut score_cutoff: f64,
    score_hint: Option<f64>,
) -> ScoreAlignment
where
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
{
    if len1 == 0 {
        return ScoreAlignment {
            score: 0.0,
            src_start: 0,
            src_end: 0,
            dest_start: 0,
            dest_end: 0,
        };
    }

    let s2_vec = s2.into_iter().collect::<Vec<_>>();

    let mut res = ScoreAlignment {
        score: 0.0,
        src_start: 0,
        src_end: len1,
        dest_start: 0,
        dest_end: len1,
    };

    let len_sum = 2 * len1;
    let mut cutoff_dist = ((len_sum as f64) * (1.0 - score_cutoff)).ceil() as usize;

    for i in 1..len1 {
        let substr_last = &s2_vec[i - 1];
        if !s1_char_set.contains(&substr_last.hash_char()) {
            continue;
        }

        let ls_ratio = comparator
            .normalized_similarity_with_args(
                s2_vec[..i].iter().cloned(),
                &indel::Args {
                    score_cutoff: WithScoreCutoff(score_cutoff),
                    score_hint,
                },
            )
            .unwrap_or(0.0);
        if ls_ratio > res.score {
            score_cutoff = ls_ratio;
            cutoff_dist = ((len_sum as f64) * (1.0 - score_cutoff)).ceil() as usize;
            res.score = ls_ratio;
            res.dest_start = 0;
            res.dest_end = i;
            if res.score == 1.0 {
                return res;
            }
        }
    }

    let window_end = len2 - len1;
    let mut scores = vec![usize::MAX; window_end + 1];
    let mut windows = vec![(0, window_end)];
    let mut new_windows = Vec::new();
    let mut best_dist = ((len_sum as f64) * (1.0 - res.score)).ceil() as usize;

    while !windows.is_empty() {
        for (start, end) in windows.drain(..) {
            if scores[start] == usize::MAX {
                let subseq = &s2_vec[start..start + len1];
                let dist = comparator.distance(subseq.iter().cloned());
                scores[start] = dist;

                if dist < best_dist {
                    best_dist = dist;
                    cutoff_dist = dist;
                    res.score = 1.0 - (dist as f64 / len_sum as f64);
                    res.dest_start = start;
                    res.dest_end = start + len1;
                    score_cutoff = res.score;
                    if dist == 0 {
                        return res;
                    }
                }
            }

            if scores[end] == usize::MAX {
                let subseq = &s2_vec[end..end + len1];
                let dist = comparator.distance(subseq.iter().cloned());
                scores[end] = dist;

                if dist < best_dist {
                    best_dist = dist;
                    cutoff_dist = dist;
                    res.score = 1.0 - (dist as f64 / len_sum as f64);
                    res.dest_start = end;
                    res.dest_end = end + len1;
                    score_cutoff = res.score;
                    if dist == 0 {
                        return res;
                    }
                }
            }

            let cell_diff = end - start;
            if cell_diff <= 1 {
                continue;
            }

            let score_start = scores[start];
            let score_end = scores[end];

            let min_val = score_start.min(score_end);
            let known_edits = score_start.abs_diff(score_end);

            // half of the cells that are not needed for known_edits can lead to a better score
            let max_score_improvement = (cell_diff.saturating_sub(known_edits / 2) / 2) * 2;
            if min_val > cutoff_dist.saturating_add(max_score_improvement) {
                continue;
            }

            let center = cell_diff / 2;
            new_windows.push((start, start + center));
            new_windows.push((start + center, end));
        }
        std::mem::swap(&mut windows, &mut new_windows);
    }

    for i in window_end + 1..len2 {
        let substr_first = &s2_vec[i];
        if !s1_char_set.contains(&substr_first.hash_char()) {
            continue;
        }

        let ls_ratio = comparator
            .normalized_similarity_with_args(
                s2_vec[i..].iter().cloned(),
                &indel::Args {
                    score_cutoff: WithScoreCutoff(score_cutoff),
                    score_hint,
                },
            )
            .unwrap_or(0.0);

        if ls_ratio > res.score {
            score_cutoff = ls_ratio;
            res.score = ls_ratio;
            res.dest_start = i;
            res.dest_end = len2;
            if res.score == 1.0 {
                return res;
            }
        }
    }

    res
}

/// `One x Many` comparisons using `partial_ratio`
pub struct PartialRatioBatchComparator<Elem1> {
    scorer: indel::BatchComparator<Elem1>,
    s1_char_set: HashSet<Hash>,
}

impl<Elem1> PartialRatioBatchComparator<Elem1>
where
    Elem1: HashableChar + Clone,
{
    pub fn new<Iter1>(s1: Iter1) -> Self
    where
        Iter1: IntoIterator<Item = Elem1>,
        Iter1::IntoIter: Clone,
    {
        let scorer = indel::BatchComparator::new(s1);
        let s1_char_set = scorer
            .scorer
            .s1
            .iter()
            .map(|c| c.hash_char())
            .collect::<HashSet<Hash>>();
        Self {
            scorer,
            s1_char_set,
        }
    }

    pub fn similarity<Iter2>(&self, s2: Iter2) -> f64
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
    {
        self.similarity_with_args(s2, &Args::default())
    }

    pub fn similarity_with_args<Iter2, CutoffType>(
        &self,
        s2: Iter2,
        args: &Args<f64, CutoffType>,
    ) -> CutoffType::Output
    where
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Elem1: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Elem1> + HashableChar + Copy,
        CutoffType: SimilarityCutoff<f64>,
    {
        let s2_iter = s2.into_iter();
        let len1 = self.scorer.scorer.s1.len();
        let len2 = s2_iter.clone().count();

        let score = if len1 <= len2 {
            partial_ratio_impl(
                &self.scorer,
                len1,
                s2_iter.clone(),
                len2,
                &self.s1_char_set,
                args.score_cutoff.cutoff().unwrap_or(0.0),
                args.score_hint,
            )
            .score
        } else {
            partial_ratio_with_args(self.scorer.scorer.s1.iter().copied(), s2_iter, args)
                .into()
                .unwrap_or(0.0)
        };

        args.score_cutoff.score(score)
    }
}

/// Computes the token ratio between two sequences.
pub fn token_ratio_with_args<Iter1, Iter2, CutoffType, CharT>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator<Item = CharT>,
    Iter2: IntoIterator<Item = CharT>,
    Iter1::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    Iter2::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    CharT: HashableChar + Clone + Ord + IsSpace + Copy,
    CutoffType: SimilarityCutoff<f64>,
{
    let score_cutoff = args.score_cutoff.cutoff().unwrap_or(0.0);

    if score_cutoff > 1.0 {
        return args.score_cutoff.score(0.0);
    }

    let tokens_a = sorted_split(s1);
    let tokens_b = sorted_split(s2);

    let decomposition = set_decomposition(tokens_a.clone(), tokens_b.clone());
    let intersect = decomposition.intersection;
    let diff_ab = decomposition.difference_ab;
    let diff_ba = decomposition.difference_ba;

    if !intersect.empty() && (diff_ab.empty() || diff_ba.empty()) {
        return args.score_cutoff.score(1.0);
    }

    let diff_ab_joined = diff_ab.join();
    let diff_ba_joined = diff_ba.join();

    let ab_len = diff_ab_joined.len();
    let ba_len = diff_ba_joined.len();
    let sect_len = intersect.length();

    let mut result = ratio_with_args(tokens_a.join(), tokens_b.join(), args)
        .into()
        .unwrap_or(0.0);

    let sect_ab_len = sect_len + usize::from(sect_len > 0) + ab_len;
    let sect_ba_len = sect_len + usize::from(sect_len > 0) + ba_len;
    let total_len = sect_ab_len + sect_ba_len;

    let cutoff_distance = score_cutoff_to_distance(score_cutoff, total_len);
    let dist_args = indel::Args::<usize, WithScoreCutoff<usize>> {
        score_cutoff: WithScoreCutoff(cutoff_distance),
        score_hint: None,
    };

    if let Some(distance) = indel::distance_with_args(diff_ab_joined, diff_ba_joined, &dist_args) {
        result = result.max(norm_distance(distance, total_len, score_cutoff));
    }

    if sect_len == 0 {
        return args.score_cutoff.score(result);
    }

    let sect_ab_dist = usize::from(sect_len > 0) + ab_len;
    let sect_ba_dist = usize::from(sect_len > 0) + ba_len;

    let sect_ab_ratio = norm_distance(sect_ab_dist, sect_len + sect_ab_len, score_cutoff);
    let sect_ba_ratio = norm_distance(sect_ba_dist, sect_len + sect_ba_len, score_cutoff);

    args.score_cutoff
        .score(result.max(sect_ab_ratio.max(sect_ba_ratio)))
}

/// Computes the partial token ratio between two sequences.
pub fn partial_token_ratio_with_args<Iter2, CutoffType, CharT>(
    s1_sorted: Vec<CharT>,
    tokens_s1: SplittedSentence<CharT>,
    s2: Iter2,
    args: &Args<f64, CutoffType>,
) -> CutoffType::Output
where
    Iter2: IntoIterator<Item = CharT>,
    Iter2::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    CharT: HashableChar + Clone + Ord + IsSpace + Copy,
    CutoffType: SimilarityCutoff<f64>,
{
    let score_cutoff = args.score_cutoff.cutoff().unwrap_or(0.0);
    if score_cutoff > 1.0 {
        return args.score_cutoff.score(0.0);
    }

    let tokens_b = sorted_split(s2);
    let decomposition = set_decomposition(tokens_s1.clone(), tokens_b.clone());

    // exit early when there is a common word in both sequences
    if !decomposition.intersection.empty() {
        return args.score_cutoff.score(1.0);
    }

    let diff_ab = decomposition.difference_ab;
    let diff_ba = decomposition.difference_ba;

    let result = partial_ratio(s1_sorted.clone(), tokens_b.join(), score_cutoff);

    // do not calculate the same partial_ratio twice
    if tokens_s1.word_count() == diff_ab.word_count()
        && tokens_b.word_count() == diff_ba.word_count()
    {
        return args.score_cutoff.score(result);
    }

    let updated_score_cutoff = score_cutoff.max(result);
    let additional_result = partial_ratio(diff_ab.join(), diff_ba.join(), updated_score_cutoff);

    args.score_cutoff.score(result.max(additional_result))
}

pub fn partial_token_ratio<Iter1, Iter2, CharT>(s1: Iter1, s2: Iter2, score_cutoff: f64) -> f64
where
    Iter1: IntoIterator<Item = CharT> + Clone,
    Iter2: IntoIterator<Item = CharT>,
    Iter1::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    Iter2::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    CharT: HashableChar + Clone + Ord + IsSpace + Copy,
{
    let tokens_s1 = sorted_split(s1.clone());
    let s1_sorted = tokens_s1.join();
    partial_token_ratio_with_args(
        s1_sorted,
        tokens_s1,
        s2,
        &Args::default().score_cutoff(score_cutoff),
    )
    .unwrap_or(0.0)
}

/// Computes the weighted ratio between two sequences.
pub fn wratio_with_args<Iter1, Iter2, CutoffType, CharT>(
    s1: Iter1,
    s2: Iter2,
    args: &Args<f64, CutoffType>,
) -> CutoffType::Output
where
    Iter1: IntoIterator<Item = CharT>,
    Iter2: IntoIterator<Item = CharT>,
    Iter1::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    Iter2::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    CharT: HashableChar + Clone + Ord + IsSpace + Copy,
    CutoffType: SimilarityCutoff<f64> + Clone + Copy,
{
    const UNBASE_SCALE: f64 = 0.95;

    let initial_cutoff = args.score_cutoff.cutoff().unwrap_or(0.0);
    if initial_cutoff > 1.0 {
        return args.score_cutoff.score(0.0);
    }

    let s1_iter = s1.into_iter();
    let s2_iter = s2.into_iter();

    let len1 = s1_iter.clone().count();
    let len2 = s2_iter.clone().count();

    // compatibility with fuzzywuzzy
    if len1 == 0 || len2 == 0 {
        return args.score_cutoff.score(0.0);
    }

    let len_ratio = if len1 > len2 {
        len1 as f64 / len2 as f64
    } else {
        len2 as f64 / len1 as f64
    };

    let mut end_ratio = ratio_with_args(s1_iter.clone(), s2_iter.clone(), args)
        .into()
        .unwrap_or(0.0);

    if len_ratio < 1.5 {
        let adjusted_cutoff = initial_cutoff.max(end_ratio) / UNBASE_SCALE;
        let new_args = args.clone().score_cutoff(adjusted_cutoff);
        let token_ratio_raw: Option<f64> =
            token_ratio_with_args(s1_iter.clone(), s2_iter.clone(), &new_args);
        let token_ratio_score = token_ratio_raw.unwrap_or(0.0);
        end_ratio = end_ratio.max(token_ratio_score * UNBASE_SCALE);
        return args.score_cutoff.score(end_ratio);
    }

    let partial_scale = if len_ratio <= 8.0 { 0.9 } else { 0.6 };

    let adjusted_cutoff = initial_cutoff.max(end_ratio) / partial_scale;
    let new_args = args.clone().score_cutoff(adjusted_cutoff);
    let partial_ratio_raw: Option<f64> =
        partial_ratio_with_args(s1_iter.clone(), s2_iter.clone(), &new_args);
    let partial_ratio_score = partial_ratio_raw.unwrap_or(0.0);
    end_ratio = end_ratio.max(partial_ratio_score * partial_scale);

    let final_cutoff = initial_cutoff.max(end_ratio) / UNBASE_SCALE;
    let new_args = args.clone().score_cutoff(final_cutoff);

    let tokens_a = sorted_split(s1_iter.clone());
    let s1_sorted = tokens_a.join();

    let partial_token_ratio_raw: Option<f64> =
        partial_token_ratio_with_args(s1_sorted, tokens_a, s2_iter, &new_args);
    let partial_token_ratio_score = partial_token_ratio_raw.unwrap_or(0.0);

    end_ratio = end_ratio.max(partial_token_ratio_score * UNBASE_SCALE * partial_scale);

    args.score_cutoff.score(end_ratio)
}

pub fn wratio<Iter1, Iter2, CharT>(s1: Iter1, s2: Iter2, score_cutoff: f64) -> f64
where
    Iter1: IntoIterator<Item = CharT>,
    Iter2: IntoIterator<Item = CharT>,
    Iter1::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    Iter2::IntoIter: Clone + DoubleEndedIterator<Item = CharT>,
    CharT: HashableChar + Clone + Ord + IsSpace + Copy,
{
    wratio_with_args(s1, s2, &Args::default().score_cutoff(score_cutoff)).unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_delta {
        ($x:expr, $y:expr) => {
            assert!(($x - $y).abs() < 0.0001, "{} != {}", $x, $y);
        };
    }

    #[test]
    fn ratio_basics() {
        assert_delta!(1.0, ratio("test".chars(), "test".chars()));
        assert_delta!(0.0, ratio("".chars(), "test".chars()));
        assert_delta!(1.0, ratio("".chars(), "".chars()));
    }

    // https://github.com/rapidfuzz/RapidFuzz/issues/206
    #[test]
    fn issue206_ratio_cutoff() {
        let str1 = "South Korea";
        let str2 = "North Korea";
        let score = ratio(str1.chars(), str2.chars());

        assert_eq!(
            None,
            ratio_with_args(
                str1.chars(),
                str2.chars(),
                &Args::default().score_cutoff(score + 0.0001)
            )
        );

        let accepted = ratio_with_args(
            str1.chars(),
            str2.chars(),
            &Args::default().score_cutoff(score - 0.0001),
        )
        .unwrap();
        assert_delta!(score, accepted);
    }

    #[test]
    fn partial_ratio_basics() {
        assert_delta!(1.0, partial_ratio("abcd".chars(), "xxabcdyy".chars(), 0.0));
        assert_delta!(1.0, partial_ratio("".chars(), "".chars(), 0.0));
        assert_delta!(0.0, partial_ratio("".chars(), "abc".chars(), 0.0));
    }

    #[test]
    fn partial_ratio_alignment_swapped() {
        let s1 = "abcd";
        let s2 = "xxabcdyy";
        let a1 = partial_ratio_alignment(
            s1.chars(),
            s1.chars().count(),
            s2.chars(),
            s2.chars().count(),
            &Args::default(),
        );
        let a2 = partial_ratio_alignment(
            s2.chars(),
            s2.chars().count(),
            s1.chars(),
            s1.chars().count(),
            &Args::default(),
        );

        assert_delta!(a1.score, a2.score);
        assert_eq!(a1.src_start, a2.dest_start);
        assert_eq!(a1.src_end, a2.dest_end);
        assert_eq!(a1.dest_start, a2.src_start);
        assert_eq!(a1.dest_end, a2.src_end);
    }

    #[test]
    fn token_ratio_basics() {
        let score = token_ratio_with_args(
            "fuzzy wuzzy was a bear".chars(),
            "wuzzy fuzzy was a bear".chars(),
            &Args::default(),
        );
        assert_delta!(1.0, score);
    }

    #[test]
    fn partial_token_ratio_basics() {
        let score = partial_token_ratio(
            "fuzzy wuzzy was a bear".chars(),
            "wuzzy fuzzy was a hare".chars(),
            0.0,
        );
        assert!(score > 0.85);
    }

    #[test]
    fn wratio_basics() {
        let score = wratio(
            "fuzzy wuzzy was a bear".chars(),
            "wuzzy fuzzy was a bear".chars(),
            0.0,
        );
        assert!(score > 0.9);
    }
}
