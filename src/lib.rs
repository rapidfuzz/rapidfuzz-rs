extern crate difflib;

pub mod utils;
pub mod fuzz;
pub mod levenshtein;
pub mod process;
mod details;

#[cfg(test)]
mod tests {
    use fuzz;
    use utils;

    #[test]
    fn test_ratio() {
        assert_eq!(fuzz::ratio("new york mets", "new york mets", 0.0), 100.0);
        assert_ne!(fuzz::ratio("new york mets", "new YORK mets", 0.0), 100.0);
    }

    #[test]
    fn test_default_process() {
        assert_eq!(utils::default_process("\0new\0YORK(mets)  "), "new york mets");
    }

    #[test]
    fn test_token_sort_ratio() {
        assert_eq!(fuzz::token_sort_ratio("new york mets", "new york mets", 0.0), 100.0);
        assert_eq!(fuzz::token_sort_ratio("york new mets", "new york mets", 0.0), 100.0);
    }

    #[test]
    fn test_token_set_ratio() {
        assert_eq!(fuzz::token_set_ratio("new york mets", "new york mets", 0.0), 100.0);
        assert_eq!(fuzz::token_set_ratio("york new mets", "new york mets", 0.0), 100.0);
        assert_eq!(fuzz::token_set_ratio("york new new mets", "new york mets", 0.0), 100.0);
    }
}