use std::cmp::Ord;
use std::fmt::Debug;
use std::vec::Vec;

use crate::details::splitted_sentence::{is_space, IsSpace, SplittedSentence};
use crate::HashableChar;

#[derive(Debug, Clone)]
pub struct DecomposedSet<CharT> {
    pub difference_ab: SplittedSentence<CharT>,
    pub difference_ba: SplittedSentence<CharT>,
    pub intersection: SplittedSentence<CharT>,
}

/// Computes the decomposition of two splitted sentences into their intersection and differences.
///
/// This function mirrors the logic of the C++ version:
/// - Dedupe both `a` and `b`
/// - Compute intersection and differences
///
/// # Parameters
/// - `a`: a `SplittedSentence<CharT>`
/// - `b`: a `SplittedSentence<CharT>`
///
/// # Returns
/// - `DecomposedSet<CharT>` containing difference_ab, difference_ba, and intersection
///
/// # Requirements
/// `CharT` must implement `IsSpace`, `HashableChar`, `Copy`, and `Ord` to ensure tokens are deduplicated and searchable.
pub fn set_decomposition<CharT>(
    mut a: SplittedSentence<CharT>,
    mut b: SplittedSentence<CharT>,
) -> DecomposedSet<CharT>
where
    CharT: IsSpace + HashableChar + Copy + Ord,
{
    // Deduplicate both splitted sentences
    a.dedupe();
    b.dedupe();

    // difference_ba initially contains all words from b
    let mut difference_ba_tokens = b.words().clone();
    let mut intersection_tokens = Vec::new();
    let mut difference_ab_tokens = Vec::new();

    // For each token in a, check if it exists in difference_ba_tokens
    for current_a in a.words() {
        if let Some(pos) = difference_ba_tokens
            .iter()
            .position(|word| word == current_a)
        {
            // Found common token, move it to intersection
            difference_ba_tokens.remove(pos);
            intersection_tokens.push(current_a.clone());
        } else {
            // Token does not exist in b, add to difference_ab
            difference_ab_tokens.push(current_a.clone());
        }
    }

    DecomposedSet {
        difference_ab: SplittedSentence::new(difference_ab_tokens),
        difference_ba: SplittedSentence::new(difference_ba_tokens),
        intersection: SplittedSentence::new(intersection_tokens),
    }
}

/// Splits an input iterator into tokens based on whitespace, sorts them, and returns a `SplittedSentence`.
///
/// # Parameters
/// - `input`: An iterator over the input sequence.
///
/// # Returns
/// - A `SplittedSentence` containing sorted tokens.
///
/// # Notes
/// - Tokens are split based on whitespace characters determined by the `is_space` function.
/// - The function collects tokens into a vector of ranges or slices, sorts them, and constructs a `SplittedSentence`.
pub fn sorted_split<Iter, CharT>(input: Iter) -> SplittedSentence<CharT>
where
    Iter: IntoIterator<Item = CharT>,
    Iter::IntoIter: Clone + Iterator<Item = CharT>,
    CharT: IsSpace + HashableChar + Copy + Ord,
{
    let mut splitted: Vec<Vec<CharT>> = Vec::new();
    let mut iter = input.into_iter().peekable();

    while let Some(&ch) = iter.peek() {
        // Skip over any whitespace characters
        if is_space(ch) {
            iter.next();
            continue;
        }

        // Collect the token
        let mut token = Vec::new();
        while let Some(&ch) = iter.peek() {
            if is_space(ch) {
                break;
            }
            token.push(ch);
            iter.next();
        }

        if !token.is_empty() {
            splitted.push(token);
        }
    }

    // Sort the tokens
    splitted.sort();

    // Construct a SplittedSentence from the sorted tokens
    SplittedSentence::new(splitted)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScoreAlignment {
    pub score: f64,
    pub src_start: usize,
    pub src_end: usize,
    pub dest_start: usize,
    pub dest_end: usize,
}

#[derive(Default, Copy, Clone)]
pub struct NoScoreCutoff;
#[derive(Default, Copy, Clone)]
pub struct WithScoreCutoff<T>(pub T);

pub trait DistanceCutoff<T>
where
    T: Copy,
{
    type Output: Copy + Into<Option<T>> + PartialEq + Debug;

    fn cutoff(&self) -> Option<T>;
    fn score(&self, raw: T) -> Self::Output;
}

impl<T> DistanceCutoff<T> for NoScoreCutoff
where
    T: Copy + PartialEq + Debug,
{
    type Output = T;

    fn cutoff(&self) -> Option<T> {
        None
    }

    fn score(&self, raw: T) -> Self::Output {
        raw
    }
}

impl<T> DistanceCutoff<T> for WithScoreCutoff<T>
where
    T: Copy + PartialOrd + Debug,
{
    type Output = Option<T>;

    fn cutoff(&self) -> Option<T> {
        Some(self.0)
    }

    fn score(&self, raw: T) -> Self::Output {
        (raw <= self.0).then_some(raw)
    }
}

pub trait SimilarityCutoff<T>
where
    T: Copy,
{
    type Output: Copy + Into<Option<T>> + PartialEq + Debug;
    type AlignmentOutput: Copy + Into<Option<ScoreAlignment>> + PartialEq + Debug;

    fn cutoff(&self) -> Option<T>;
    fn score(&self, raw: T) -> Self::Output;
    fn alignment(&self, raw: Option<ScoreAlignment>) -> Self::AlignmentOutput;
}

impl<T> SimilarityCutoff<T> for NoScoreCutoff
where
    T: Copy + PartialEq + Debug,
{
    type Output = T;
    type AlignmentOutput = ScoreAlignment;

    fn cutoff(&self) -> Option<T> {
        None
    }

    fn score(&self, raw: T) -> Self::Output {
        raw
    }

    fn alignment(&self, raw: Option<ScoreAlignment>) -> Self::AlignmentOutput {
        raw.expect("alignment must be present when no score cutoff is set")
    }
}

impl<T> SimilarityCutoff<T> for WithScoreCutoff<T>
where
    T: Copy + PartialOrd + Debug,
{
    type Output = Option<T>;
    type AlignmentOutput = Option<ScoreAlignment>;

    fn cutoff(&self) -> Option<T> {
        Some(self.0)
    }

    fn score(&self, raw: T) -> Self::Output {
        (raw >= self.0).then_some(raw)
    }

    fn alignment(&self, raw: Option<ScoreAlignment>) -> Self::AlignmentOutput {
        raw
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_decomposition() {
        let s1_tokens = vec![
            vec!['f', 'u', 'z', 'z', 'y'],
            vec!['w', 'u', 'z', 'z', 'y'],
            vec!['w', 'a', 's'],
        ];
        let s2_tokens = vec![
            vec!['f', 'u', 'z', 'z', 'y'],
            vec!['f', 'u', 'z', 'z', 'y'],
            vec!['b', 'e', 'a', 'r'],
        ];
        let s1 = SplittedSentence::new(s1_tokens);
        let s2 = SplittedSentence::new(s2_tokens);

        let result = set_decomposition(s1, s2);

        // After dedupe:
        // s1 words: fuzzy, wuzzy, was
        // s2 words: fuzzy, bear
        // intersection: fuzzy
        // difference_ab: wuzzy, was
        // difference_ba: bear

        assert_eq!(
            result.intersection.words(),
            &vec![vec!['f', 'u', 'z', 'z', 'y']]
        );
        assert_eq!(
            result.difference_ab.words(),
            &vec![vec!['w', 'u', 'z', 'z', 'y'], vec!['w', 'a', 's']]
        );
        assert_eq!(
            result.difference_ba.words(),
            &vec![vec!['b', 'e', 'a', 'r']]
        );
    }
}
