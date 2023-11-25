use std::ops::{BitAnd, Shl};

pub struct BitMatrix<T> {
    rows: usize,
    cols: usize,
    matrix: Vec<T>,
}

impl<T> BitMatrix<T>
where
    T: Clone,
{
    pub fn new(rows: usize, cols: usize, val: T) -> Self {
        Self {
            rows,
            cols,
            matrix: vec![val; rows * cols],
        }
    }

    #[allow(dead_code)]
    pub const fn rows(&self) -> usize {
        self.rows
    }

    #[allow(dead_code)]
    pub const fn cols(&self) -> usize {
        self.cols
    }

    pub fn get(&self, row: usize, col: usize) -> &T {
        debug_assert!(row < self.rows);
        debug_assert!(col < self.cols);
        &self.matrix[row * self.cols + col]
    }

    pub fn get_mut(&mut self, row: usize, col: usize) -> &mut T {
        debug_assert!(row < self.rows);
        debug_assert!(col < self.cols);
        &mut self.matrix[row * self.cols + col]
    }
}

pub struct ShiftedBitMatrix<T> {
    matrix: BitMatrix<T>,
    offsets: Vec<isize>,
}

impl<T> ShiftedBitMatrix<T>
where
    T: Copy + From<u8> + Shl<usize, Output = T> + BitAnd<T, Output = T> + PartialEq<T>,
{
    pub fn new(rows: usize, cols: usize, val: T) -> Self {
        Self {
            matrix: BitMatrix::<T>::new(rows, cols, val),
            offsets: vec![0; rows],
        }
    }

    #[allow(dead_code)]
    pub fn test_bit(&self, row: usize, mut col: usize, default: bool) -> bool {
        let offset = self.offsets[row];

        if offset < 0 {
            col += (-offset) as usize;
        } else if col >= offset as usize {
            col -= offset as usize;
        }
        // bit on the left of the band
        else {
            return default;
        }

        let word_size = std::mem::size_of::<T>() * 8;
        let col_word = col / word_size;
        let col_mask = T::from(1) << (col % word_size);

        (*self.matrix.get(row, col_word) & col_mask) != T::from(0)
    }

    #[allow(dead_code)]
    pub fn get(&self, row: usize, col: usize) -> &T {
        self.matrix.get(row, col)
    }

    pub fn get_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.matrix.get_mut(row, col)
    }

    pub fn set_offset(&mut self, row: usize, offset: isize) {
        self.offsets[row] = offset;
    }
}

impl<T> Default for ShiftedBitMatrix<T>
where
    T: Copy + From<u8> + Shl<usize, Output = T> + BitAnd<T, Output = T> + PartialEq<T>,
{
    fn default() -> Self {
        Self::new(0, 0, T::from(0))
    }
}
