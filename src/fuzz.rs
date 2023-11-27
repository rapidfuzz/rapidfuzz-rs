use crate::distance::indel;
use crate::HashableChar;

/// Returns a simple ratio between two strings or `None` if `ratio < score_cutoff`
///
/// # Example
/// ```
/// use rapidfuzz::fuzz;
/// /// score is Some(0.9655)
/// let score = fuzz::ratio("this is a test".chars(), "this is a test!".chars(), None, None);
/// ```
///
pub fn ratio<Iter1, Iter2, Elem1, Elem2, ScoreCutoff, ScoreHint>(
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
    indel::normalized_similarity(s1, s2, score_cutoff, score_hint)
}

/// `One x Many` comparisons using `ratio`
///
/// # Examples
///
/// ```
/// use rapidfuzz::fuzz;
///
/// let scorer = fuzz::RatioBatchComparator::new("this is a test".chars());
/// /// score is Some(0.9655)
/// let score = scorer.similarity("this is a test!".chars(), None, None);
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
    pub fn similarity<Iter2, Elem2, ScoreCutoff, ScoreHint>(
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
        self.scorer
            .normalized_similarity(s2, score_cutoff, score_hint)
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
        assert_delta!(Some(1.0), ratio(S1.chars(), S1.chars(), None, None));
        assert_delta!(Some(1.0), ratio("test".chars(), "test".chars(), None, None));
        assert_delta!(Some(1.0), ratio(S8.chars(), S8.chars(), None, None));
        assert_delta!(Some(1.0), ratio(S9.chars(), S9.chars(), None, None));
    }

    #[test]
    fn test_partial_ratio() {
        //assert_delta!(Some(1.0), partial_ratio(S1.chars(), S1.chars(), None, None));
        assert_delta!(Some(0.65), ratio(S1.chars(), S3.chars(), None, None));
        //assert_delta!(Some(1.0), partial_ratio(S1.chars(), S3.chars(), None, None));
    }

    #[test]
    fn two_empty_strings() {
        assert_delta!(Some(1.0), ratio("".chars(), "".chars(), None, None));
    }

    #[test]
    fn first_string_empty() {
        assert_delta!(Some(0.0), ratio("test".chars(), "".chars(), None, None));
    }

    #[test]
    fn second_string_empty() {
        assert_delta!(Some(0.0), ratio("".chars(), "test".chars(), None, None));
    }

    // https://github.com/maxbachmann/RapidFuzz/issues/206
    #[test]
    fn issue206() {
        let str1 = "South Korea";
        let str2 = "North Korea";

        {
            let score = ratio(str1.chars(), str2.chars(), None, None).expect("no score_cutoff");

            assert_eq!(
                None,
                ratio(str1.chars(), str2.chars(), Some(score + 0.0001), None)
            );
            assert_delta!(
                Some(score),
                ratio(str1.chars(), str2.chars(), Some(score - 0.0001), None)
            );
        }
    }

    // https://github.com/maxbachmann/RapidFuzz/issues/210
    #[test]
    fn issue210() {
        let str1 = "bc";
        let str2 = "bca";

        {
            let score = ratio(str1.chars(), str2.chars(), None, None).expect("no score_cutoff");

            assert_eq!(
                None,
                ratio(str1.chars(), str2.chars(), Some(score + 0.0001), None)
            );
            assert_delta!(
                Some(score),
                ratio(str1.chars(), str2.chars(), Some(score - 0.0001), None)
            );
        }
    }
}
