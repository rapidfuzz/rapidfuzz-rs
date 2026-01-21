/**
Tuple like object describing the position of the compared strings in
src and dest.

It indicates that the score has been calculated between
src[src_start:src_end] and dest[dest_start:dest_end]
*/
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct ScoreAlignment {
    pub score: f64,
    pub src_start: usize,
    pub src_end: usize,
    pub dest_start: usize,
    pub dest_end: usize,
}
