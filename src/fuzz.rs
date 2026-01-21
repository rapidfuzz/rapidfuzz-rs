use crate::common::{NoScoreCutoff, SimilarityCutoff, WithScoreCutoff};
use crate::details::distance::MetricUsize;
use crate::distance::common::ScoreAlignment;
use crate::distance::indel;
use crate::{Hash, HashableChar};
use std::cmp;
use std::collections::HashSet;

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

/// Returns a simple ratio between two strings or `None` if `ratio < score_cutoff`
///
/// # Example
/// ```
/// use rapidfuzz::fuzz;
/// /// score is 0.9655
/// let score = fuzz::ratio("this is a test".chars(), "this is a test!".chars());
/// ```
///
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
///
/// # Examples
///
/// ```
/// use rapidfuzz::fuzz;
///
/// let scorer = fuzz::RatioBatchComparator::new("this is a test".chars());
/// /// score is 0.9655
/// let score = scorer.similarity("this is a test!".chars());
/// ```
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

    /// Similarity calculated similar to [`ratio`]
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

/// Searches for the optimal alignment of the shorter string in the
/// longer string and returns the fuzz.ratio for this alignment.
///
/// # Example
/// ```
/// use rapidfuzz::fuzz;
/// /// score is 1.0
/// let score = fuzz::partial_ratio("this is a test".chars(), "this is a test!".chars());
/// ```
///
pub fn partial_ratio<Iter1, Iter2>(s1: Iter1, s2: Iter2) -> f64
where
    Iter1: IntoIterator,
    Iter1::IntoIter: DoubleEndedIterator + Clone,
    Iter2: IntoIterator,
    Iter2::IntoIter: DoubleEndedIterator + Clone,
    Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
    Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
{
    partial_ratio_with_args(s1, s2, &Args::default())
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

    let score = alignment.into().map_or(0.0, |alignment| alignment.score);
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

    if len1 > len2 {
        let res = partial_ratio_alignment(s2_iter.clone(), len2, s1_iter.clone(), len1, args);
        return args.score_cutoff.alignment(res.into().map(|mut alignment| {
            std::mem::swap(&mut alignment.src_start, &mut alignment.dest_start);
            std::mem::swap(&mut alignment.src_end, &mut alignment.dest_end);
            alignment
        }));
    }

    if len1 == 0 && len2 == 0 {
        return args.score_cutoff.alignment(Some(ScoreAlignment {
            score: 1.0,
            src_start: 0,
            src_end: 0,
            dest_start: 0,
            dest_end: 0,
        }));
    }

    let s1_char_set = s1_iter
        .clone()
        .map(|c| c.hash_char())
        .collect::<HashSet<_>>();

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

/**
implementation of partial_ratio for needles <= 64. assumes len(s1) <= len(s2)
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
    let mut cutoff_dist = (len_sum as f64 * (1.0 - score_cutoff)).ceil() as usize;

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
            cutoff_dist = (len_sum as f64 * (1.0 - score_cutoff)).ceil() as usize;
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
    let mut best_dist = (len_sum as f64 * (1.0 - res.score)).ceil() as usize;

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

            let min_val = cmp::min(score_start, score_end);
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
///
/// # Examples
///
/// ```
/// use rapidfuzz::fuzz;
///
/// let scorer = fuzz::PartialRatioBatchComparator::new("this is a test".chars());
/// /// score is 1.0
/// let score = scorer.similarity("this is a test!".chars());
/// ```
pub struct PartialRatioBatchComparator<Elem1> {
    scorer: indel::BatchComparator<Elem1>,
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
        Self {
            scorer: indel::BatchComparator::new(s1),
        }
    }

    /// Similarity calculated similar to [`partial_ratio`]
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

        let res = partial_ratio_impl(
            &self.scorer,
            self.scorer.scorer.s1.len(),
            s2_iter.clone(),
            s2_iter.count(),
            &self
                .scorer
                .scorer
                .s1
                .iter()
                .map(|c| c.hash_char())
                .collect(),
            args.score_cutoff.cutoff().unwrap_or(0.0),
            args.score_hint,
        );

        let alignment = args.score_cutoff.alignment(Some(res));
        let score = alignment.into().map_or(0.0, |alignment| alignment.score);
        args.score_cutoff.score(score)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static S1: &str = "new york mets";
    static S3: &str = "the wonderful new york mets";
    //static S4: &str = "new york mets vs atlanta braves";
    //static S5: &str = "atlanta braves vs new york mets";
    //static S7: &str = "new york city mets - atlanta braves";
    // test silly corner cases
    static S8: &str = "{";
    static S9: &str = "{a";
    //static S10: &str = "a{";
    //static S10A: &str = "{b";

    macro_rules! assert_delta {
        ($x:expr, $y:expr) => {
            match ($x, $y) {
                (None, None) => {}
                (Some(val1), Some(val2)) => {
                    if (val1 - val2).abs() > 0.0001 {
                        panic!("{:?} != {:?}", $x, $y);
                    }
                }
                (_, _) => panic!("{:?} != {:?}", $x, $y),
            }
        };
    }

    #[test]
    fn test_equal() {
        assert_delta!(
            Some(1.0),
            Some(ratio_with_args(S1.chars(), S1.chars(), &Args::default()))
        );
        assert_delta!(
            Some(1.0),
            Some(ratio_with_args(
                "test".chars(),
                "test".chars(),
                &Args::default()
            ))
        );
        assert_delta!(
            Some(1.0),
            Some(ratio_with_args(S8.chars(), S8.chars(), &Args::default()))
        );
        assert_delta!(
            Some(1.0),
            Some(ratio_with_args(S9.chars(), S9.chars(), &Args::default()))
        );
    }

    #[test]
    fn test_partially_equal() {
        //assert_delta!(Some(1.0), partial_ratio(S1.chars(), S1.chars(), None, None));
        assert_delta!(
            Some(0.65),
            Some(ratio_with_args(S1.chars(), S3.chars(), &Args::default()))
        );
        //assert_delta!(Some(1.0), partial_ratio(S1.chars(), S3.chars(), None, None));
    }

    #[test]
    fn two_empty_strings() {
        assert_delta!(
            Some(1.0),
            Some(ratio_with_args("".chars(), "".chars(), &Args::default()))
        );

        assert_delta!(
            Some(1.0),
            Some(partial_ratio_with_args(
                "".chars(),
                "".chars(),
                &Args::default()
            ))
        );
    }

    #[test]
    fn first_string_empty() {
        assert_delta!(
            Some(0.0),
            Some(ratio_with_args(
                "test".chars(),
                "".chars(),
                &Args::default()
            ))
        );
    }

    #[test]
    fn second_string_empty() {
        assert_delta!(
            Some(0.0),
            Some(ratio_with_args(
                "".chars(),
                "test".chars(),
                &Args::default()
            ))
        );
    }

    // https://github.com/rapidfuzz/RapidFuzz/issues/206
    #[test]
    fn issue206() {
        let str1 = "South Korea";
        let str2 = "North Korea";

        {
            let score = ratio(str1.chars(), str2.chars());

            assert_eq!(
                None,
                ratio_with_args(
                    str1.chars(),
                    str2.chars(),
                    &Args::default().score_cutoff(score + 0.0001)
                )
            );
            assert_delta!(
                Some(score),
                ratio_with_args(
                    str1.chars(),
                    str2.chars(),
                    &Args::default().score_cutoff(score - 0.0001)
                )
            );
        }

        {
            let score = partial_ratio(str1.chars(), str2.chars());

            assert_eq!(
                None,
                partial_ratio_with_args(
                    str1.chars(),
                    str2.chars(),
                    &Args::default().score_cutoff(score + 0.0001)
                )
            );
            assert_delta!(
                Some(score),
                partial_ratio_with_args(
                    str1.chars(),
                    str2.chars(),
                    &Args::default().score_cutoff(score - 0.0001)
                )
            );
        }
    }

    // https://github.com/rapidfuzz/RapidFuzz/issues/210
    #[test]
    fn issue210() {
        let str1 = "bc";
        let str2 = "bca";

        {
            let score = ratio(str1.chars(), str2.chars());

            assert_eq!(
                None,
                ratio_with_args(
                    str1.chars(),
                    str2.chars(),
                    &Args::default().score_cutoff(score + 0.0001)
                )
            );
            assert_delta!(
                Some(score),
                ratio_with_args(
                    str1.chars(),
                    str2.chars(),
                    &Args::default().score_cutoff(score - 0.0001)
                )
            );
        }
    }

    #[test]
    fn test_partial_ratio() {
        let s1 = "this is a test";
        let s2 = "this is a test!";
        let result = partial_ratio(s1.chars(), s2.chars());

        assert_delta!(Some(1.0), Some(result));
    }

    #[test]
    fn test_partial_ratio_issue76() {
        let s1 = "physics 2 vid";
        let s2 = "study physics physics 2 video";

        let result = partial_ratio(s1.chars(), s2.chars());

        assert_delta!(Some(1.0), Some(result))
    }

    #[test]
    fn test_partial_ratio_issue90() {
        let s1 = "ax b";
        let s2 = "a b a c b";

        let result = partial_ratio(s1.chars(), s2.chars());

        assert_delta!(Some(0.8571428), Some(result))
    }

    #[test]
    fn test_partial_ratio_issue138() {
        let s1 = &"a".repeat(65);
        let s2 = &format!("a{}{}", char::from_u32(256).unwrap(), "a".repeat(63));

        let result = partial_ratio(s1.chars(), s2.chars());

        assert_delta!(Some(0.9922481), Some(result));
    }

    #[test]
    fn test_partial_ratio_issue219() {
        let s1 = "\
            TTAGCGCTACCGGTCGCCACCATGGTTTTCTAAGGGGAGGCCGTCATCAAAAGAGTTCATGTAGCACGAAGTCCACCTTTGAAGGATCGATGAATG\
            GCCATGAATTCGAAATCGAGGGGAGGGCGAGAGAGGGCCGGCCTTACGAGGGCACACCCAAACTGCCAAACTGAAAGTGACCAAAGGCGGCCCGTT\
            ACCATTCTCCTGGGACATACTGTAAGTGCATGGCACCACGCTCTATTTCTTAAAAAAAGTGTAGGGTCTGGCGCCCTCGGGGGCGGCTTAGGAAAA\
            GAGGCCTGACCAATTTTTGTCTCTTATAGGTCACCACAGTTCATGTACGGAAGCAGAGCGTTCACGAAGCACCCAGCTGACATCCCGGACTACTAT\
            GACAGAGCTTCCCGGAAGGACTCAAGTGGGAGCGGGTCATGAACTTCGAGGACGGTGGGGCAGTGACTGTGACACAGGACACCAGCCTGAAGATGG\
            AACTCTTATCTACAAAGTAAAGCTAAGAGGAACCAACTTCCCGCCAGATGGGCCCGTTATGCAAAAGAAAACGATGGGGTGGGAAGCTTCTGCAGA\
            GCGCCTTTACCCCGAGGATGGCGTCCTTAAGGGGGATATCAAAATGGCGCTACGCCTTAAGGATGGAGGCAGATATTTGGCAGACTTCAAAACAAC\
            ATTACAAGGCGAAGAAGCCAGTCCAGATGCCTGGAGCTTGCAATGGTAAGCACCTCTGCCTGCCCCGCTAGTTGGGTGTGAGTGGCCCAGGCAGCC\
            GCCTGCATTTAGCTCTAGCCGGGGTACGGGTGCCCCTTGATGCCTGAGGCCTCTCCTGTGGCTGAGGCGACTGGCCCAGAGTCTGGGTCTCCTCGA\
            GGGTGGCCATCTGGCGTCACCTGTCATCTGCCACCTCTGACCCCTGCCTCTCTCCTCACAGTTGACCGGAAGCTCGACATAACGAGTCACAACGAG\
            GACTACACAGTTGTCGAGCAGTACGAACGTTCCGAGGGTCGACACTCAACTGGCAGGATGGATGAGCTTTTACAAAGGGCGGGGGCGGAGGAAGCG\
            GAGGAGGAGGAAGTGGTGGAGGAGGCTCGAAAGGTAAGTATCAGGGTTGCAGCGTTTCTCTGACCTCATATTCCAATGGATGTGTGAGAAGCATAG\
            TGAGATCCGTTTACCCCTTTTGCTCAATTCTCACGTGGCTGTAGTCGTGTTTATAAGTCTGATCGTAATGGCAGCTTGGTCTGCGTGCCTTGAAAT\
            TGTGGCCCCCACATGCATAATAAACGATCCTCTAGCACTACTTTCTGTCGAGCCACCTCAGCGCCCGTACAGTAATGTCTACAGCGCGTCTAACCC\
            GACAAATGCGTTTCTTTCTCTCCTAGAACGAAAGATTACGGATCACAGAAACGTCTCGGAAAGTCCAAATAGAAAGAACGAGAAAGAAGAAAGTGA\
            AGGATCACAAGAGCAACTCGAAAGAAAGAGACATAAGAAGGAACTCAGAAAAGGATGACAAGTATAAAAACAAAGTGAAGAAAAGAGCGAAGAGCA\
            GAGTAGAAGCAAGAGTAAAGAGAAGAAGAGCAAATCGAAGGAAAGGTAAGTGGCTTTCAAGAACATTGGTAAAACGTCATGTGTATTGCGGTTCCA\
            TGCTTACACAAATTCGTTCGCTTGTTTTCAGGGACTCGAAACACAACAGAAACGAAGAGAAGAGAATGAGAAGCAGAAGCAAAGGAAGAGACCATG\
            AAAATGTCAAGGAAAAAGAAAAAACAGTCCGATAGCAAAGGCAAAGACCAGGAGCGGTCTCGGTCGAAGGAAAAATCTAAACAACTTGAATCAAAA\
            TCTAACGAGCATGGTAAGTTCGCGAGACACTAAGTTGATTCTTAGTGTTTAGACGTGAAACTCCCTTGGAAGGTTTAACGAATACTGTTAATATTT\
            TCAGATCACTCAAAATCCAAAAGAACCGACGGGCACAATCCCGGAGCCGTGAATGTGATATAACCAAGGAAGCACAGTTGCAATTCGAGAACAAGA\
            GAAAGAAGCAGAAGTAGAGAGATCGCTCGAGAAGAGTGAGAAGCAGAACACATGATAGAGACAGAAGCCGGTCGAAAGAATACCACCGCTACAGAG\
            AACAAGGTAAGCATGACTACTTGAGTGTAAATACGTTGTGATAGAGATGAAAAACAAAACCGAACATTACTTTGGGTAATAATTAACTTTTTTTTA\
            ATAGAATATCGGGAGAAAGGAAGGTCGAGAAGCAGAGAAAGAAGGACGCCTCAGGAAGAAGCCGTTCGAAAGACAGAAGGAGAAGGAGAAGAGATT\
            CGAAAGTTCAGAGCGTGAAGAGTCTCAATCGCGTAATAAAGACAAGTACGGGAACCAAGAAAGTAAAAGTTCCCACAGGAAGAACTCTGAAGAGCG\
            AGAAAAGTAAAAAAGGGTTTCCTGTTTTTTGCCTATTTTGGGTAAAGGGGTTGATGGAGAAACAGGTGTGTGGACTGCTGAGGAGTGAGTTAGAAT\
            AAATGGTGGTATCACTTCTTCAATGCTACTACAATGGAACAACAGTCGTTACCTGTTTTAAGTTCGTGGCGTCTTATGCTCCGGACAGGGACAGAT\
            AGGCGGTTGACAGAGAGTTAAGATCTAGTACACTGGGTTTCCTAAATGTAAGAATTGGCCCGAATCCGGCCTAATATGCGAACTTTGTGCTACCAA\
            GCGAGCGGGAAGCTAAGGGTGGGGAATTGCGGGTTTAATGGACCATCTCATGAGTCTAGCAGTTAATGTATCCTATCTTCCAAACAGGAATGTATT\
            CGAAAGAGTAGAGACCATAATTCGTCTAACAACTCAAGGAAAAGAAGGCGGAGTAGAGCCGATTCCGAACCCTTTGCTAGGACTAGATAGCACGTG\
            AACCTAGACTGTCTCTGAGACTGCGCCATTACGTCTCGATCAGTAACGATTGCATCGCGAGGCTGTGGATGTAAAACCTCTGCTGACCTTGACTGA\
            CTGAGATACAATGCCTTCAGCAATGCGTGGCAG";
        let s2 = "\
            GTAAGGGTTTCCTGTTTTTTGCCTATTTTGGGTAAAGGGGGGTTGATGGAGAAACAGGTGTGTGGACTGCTGAGGAGTGAGTTAGAATAAATGGTG\
            GTATCACTTCTTCAATGCTACAATGGAACAACAGTCGTTACCTGTTTTAAGTTCGTGGCGTCTTATGCTCCGGACAGGGACAGATAGGCGGTTAGA\
            CAGAGAGTTAAGATCTAGTACACTGGGTTTCCTAAATGTAAAAATTGGCCCGAATCCGGCCTAATATGCGAACTTTGTGCTACCAAGCGAGCGGGA\
            AGCTAAGGGTGGGGAGTGCGGGTTTAATGGACCATCTCGCAGGTCTAGCAGTTAATGTATCCTATCTTCCAAACAG";

        let expected = 0.975274725;

        let score = partial_ratio(s1.chars(), s2.chars());
        assert_delta!(Some(expected), Some(score));

        let score = partial_ratio(s2.chars(), s1.chars());
        assert_delta!(Some(expected), Some(score));
    }

    #[test]
    fn test_partial_ratio_issue257() {
        let s1 = "aaaaaaaaaaaaaaaaaaaaaaaabacaaaaaaaabaaabaaaaaaaababbbbbbbbbbabbcb";
        let s2 = "aaaaaaaaaaaaaaaaaaaaaaaababaaaaaaaabaaabaaaaaaaababbbbbbbbbbabbcb";

        let expected = 0.9846153846153847;

        let score = partial_ratio(s1.chars(), s2.chars());
        assert_delta!(Some(expected), Some(score));

        let score = partial_ratio(s2.chars(), s1.chars());
        assert_delta!(Some(expected), Some(score));
    }

    #[test]
    fn test_partial_ratio_alignment() {
        let s1 = "a certain string";
        let s2 = "certain";

        let result1 = partial_ratio_alignment(
            s2.chars(),
            s2.chars().count(),
            s1.chars(),
            s1.chars().count(),
            &Args::default(),
        );

        assert_eq!(
            result1,
            ScoreAlignment {
                score: 1.0,
                src_start: 0,
                src_end: s2.len(),
                dest_start: 2,
                dest_end: 2 + s2.len()
            }
        );

        let result2 = partial_ratio_alignment(
            s1.chars(),
            s1.chars().count(),
            s2.chars(),
            s2.chars().count(),
            &Args::default(),
        );

        assert_eq!(
            result2,
            ScoreAlignment {
                score: 1.0,
                src_start: 2,
                src_end: 2 + s2.len(),
                dest_start: 0,
                dest_end: s2.len()
            }
        );

        // assert_eq!(
        //     None,
        //     Some(partial_ratio_alignment(
        //         None,
        //         0,
        //         "test".chars(),
        //         "test".chars().count(),
        //         &Args::default()
        //     ))
        // );

        assert_eq!(
            None,
            partial_ratio_alignment(
                "test".chars(),
                "test".chars().count(),
                "tesx".chars(),
                "tesx".chars().count(),
                &Args::default().score_cutoff(0.9)
            )
        );
    }

    #[test]
    fn test_partial_ratio_alignment_issue231() {
        let str1 = "er merkantilismus förderte handle und verkehr mit teils marktkonformen, teils dirigistischen maßnahmen.";
        let str2 = "ils marktkonformen, teils dirigistischen maßnahmen. an der schwelle zum 19. jahrhundert entstand ein neu";

        let alignment = partial_ratio_alignment(
            str1.chars(),
            str1.chars().count(),
            str2.chars(),
            str2.chars().count(),
            &Args::default(),
        );

        dbg!(&alignment);

        assert!(
            (alignment.score - 0.662337662).abs() < 1e-5,
            "Expected 0.662337662, got {}",
            alignment.score
        );
        assert_eq!(alignment.src_start, 0);
        assert_eq!(alignment.src_end, 103);
        assert_eq!(alignment.dest_start, 0);
        assert_eq!(alignment.dest_end, 51);
    }

    #[test]
    fn test_partial_ratio_impl_identical() {
        let s1 = "abcd";
        let s2 = "abcd";

        let indel_comp = indel::BatchComparator::new(s1.chars());

        let result = partial_ratio_impl(
            &indel_comp,
            s1.chars().count(),
            s2.chars(),
            s2.chars().count(),
            &s1.chars().map(|c| c.hash_char()).collect(),
            0.0,
            None,
        );

        assert_eq!(result.score, 1.0);
        assert_eq!(result.src_start, 0);
        assert_eq!(result.src_end, 4);
        assert_eq!(result.dest_start, 0);
        assert_eq!(result.dest_end, 4);
    }

    #[test]
    fn test_partial_ratio_impl_substring() {
        let s1 = "bcd";
        let s2 = "abcde";

        let indel_comp = indel::BatchComparator::new(s1.chars());

        let result = partial_ratio_impl(
            &indel_comp,
            s1.chars().count(),
            s2.chars(),
            s2.chars().count(),
            &s1.chars().map(|c| c.hash_char()).collect(),
            0.0,
            None,
        );

        assert_eq!(result.score, 1.0);
        assert_eq!(result.src_start, 0);
        assert_eq!(result.src_end, 3);
        assert_eq!(result.dest_start, 1);
        assert_eq!(result.dest_end, 4);
    }
}
