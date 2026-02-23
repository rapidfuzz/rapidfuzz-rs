use crate::HashableChar;

// src/details/splitted_sentence.rs

/// Trait to determine if a character is a whitespace and to provide a space character.
pub trait IsSpace: Sized + Copy {
    /// Determines if the character is a whitespace character.
    fn is_space(&self) -> bool;

    /// Returns a space character of the same type.
    fn space() -> Self;
}

impl IsSpace for char {
    fn is_space(&self) -> bool {
        matches!(
            *self,
            '\u{0009}' // TAB
                | '\u{000A}' // LF
                | '\u{000B}' // VT
                | '\u{000C}' // FF
                | '\u{000D}' // CR
                | '\u{001C}'
                | '\u{001D}'
                | '\u{001E}'
                | '\u{001F}'
                | '\u{0020}' // SPACE
                | '\u{0085}'
                | '\u{00A0}'
                | '\u{1680}'
                | '\u{2000}'
                | '\u{2001}'
                | '\u{2002}'
                | '\u{2003}'
                | '\u{2004}'
                | '\u{2005}'
                | '\u{2006}'
                | '\u{2007}'
                | '\u{2008}'
                | '\u{2009}'
                | '\u{200A}'
                | '\u{2028}'
                | '\u{2029}'
                | '\u{202F}'
                | '\u{205F}'
                | '\u{3000}'
        )
    }

    fn space() -> Self {
        ' '
    }
}

impl IsSpace for u8 {
    fn is_space(&self) -> bool {
        matches!(
            *self,
            0x09 | 0x0A | 0x0B | 0x0C | 0x0D | 0x1C | 0x1D | 0x1E | 0x1F | 0x20
        )
    }

    fn space() -> Self {
        0x20 // ASCII space
    }
}

/// Determines if a character is considered a whitespace character.
///
/// This function now operates on any type that implements the `IsSpace` trait.
pub fn is_space<CharT: IsSpace>(ch: CharT) -> bool {
    ch.is_space()
}

/// A view into a splitted sentence, containing sorted tokens.
#[derive(Debug, Clone)]
pub struct SplittedSentence<CharT> {
    tokens: Vec<Vec<CharT>>,
}

impl<CharT> SplittedSentence<CharT>
where
    CharT: IsSpace + HashableChar + Copy + Ord,
{
    /// Creates a new `SplittedSentence` from a vector of token vectors.
    pub fn new(tokens: Vec<Vec<CharT>>) -> Self {
        SplittedSentence { tokens }
    }

    /// Removes duplicate tokens, keeping only unique tokens.
    ///
    /// Returns the number of duplicates removed.
    pub fn dedupe(&mut self) -> usize {
        let old_word_count = self.word_count();
        self.tokens.dedup(); // Removes consecutive duplicates while preserving order.
        old_word_count - self.word_count()
    }

    /// Returns the total size (number of characters plus spaces) of the splitted sentence.
    pub fn size(&self) -> usize {
        if self.tokens.is_empty() {
            return 0;
        }

        // There is a space between each word
        let mut result = self.tokens.len() - 1;
        for token in &self.tokens {
            result += token.len();
        }

        result
    }

    /// Returns the length of the splitted sentence.
    ///
    /// This is an alias for `size`.
    pub fn length(&self) -> usize {
        self.size()
    }

    /// Checks if the splitted sentence is empty.
    pub fn empty(&self) -> bool {
        self.tokens.is_empty()
    }

    /// Returns the number of words (tokens) in the splitted sentence.
    pub fn word_count(&self) -> usize {
        self.tokens.len()
    }

    /// Joins the tokens back into a single vector of characters, separated by spaces.
    pub fn join(&self) -> Vec<CharT> {
        if self.tokens.is_empty() {
            return Vec::new();
        }

        let mut joined = Vec::with_capacity(self.size());
        joined.extend(&self.tokens[0]);

        for token in self.tokens.iter().skip(1) {
            joined.push(CharT::space());
            joined.extend(token);
        }

        joined
    }

    /// Returns a reference to the internal tokens.
    pub fn words(&self) -> &Vec<Vec<CharT>> {
        &self.tokens
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_splitted_sentence_char() {
        let tokens = vec![
            vec!['f', 'u', 'z', 'z', 'y'],
            vec!['w', 'u', 'z', 'z', 'y'],
            vec!['w', 'a', 's'],
            vec!['a'],
            vec!['b', 'e', 'a', 'r'],
        ];
        let mut splitted = SplittedSentence::new(tokens.clone());
        // 'fuzzy wuzzy was a bear' has 5 + 1 + 5 + 1 + 3 + 1 + 1 + 1 + 4 = 22 characters
        assert_eq!(splitted.size(), 22);

        let removed = splitted.dedupe();
        // All tokens are unique, so dedupe should remove 0
        assert_eq!(removed, 0);
        assert_eq!(splitted.word_count(), 5);

        let joined = splitted.join();
        assert_eq!(
            joined,
            vec![
                'f', 'u', 'z', 'z', 'y', ' ', 'w', 'u', 'z', 'z', 'y', ' ', 'w', 'a', 's', ' ',
                'a', ' ', 'b', 'e', 'a', 'r'
            ]
        );
    }

    #[test]
    fn test_splitted_sentence_u8() {
        let tokens = vec![
            vec![102, 117, 122, 122, 121], // "fuzzy"
            vec![119, 117, 122, 122, 121], // "wuzzy"
            vec![119, 97, 115],            // "was"
            vec![97],                      // "a"
            vec![98, 101, 97, 114],        // "bear"
        ];
        let mut splitted = SplittedSentence::new(tokens.clone());
        // 'fuzzy wuzzy was a bear' has 5 + 1 + 5 + 1 + 3 + 1 + 1 + 1 + 4 = 22 characters
        assert_eq!(splitted.size(), 22);

        let removed = splitted.dedupe();
        // All tokens are unique, so dedupe should remove 0
        assert_eq!(removed, 0);
        assert_eq!(splitted.word_count(), 5);

        let joined = splitted.join();
        assert_eq!(
            joined,
            vec![
                102, 117, 122, 122, 121, 32, // "fuzzy "
                119, 117, 122, 122, 121, 32, // "wuzzy "
                119, 97, 115, 32, // "was "
                97, 32, // "a "
                98, 101, 97, 114 // "bear"
            ]
        );
    }
}
