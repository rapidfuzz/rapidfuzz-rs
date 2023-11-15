pub(crate) fn norm_sim_to_norm_dist(score_cutoff: f64) -> f64 {
    let imprecision = 0.00001;
    (1.0 - score_cutoff + imprecision).min(1.0)
}

pub enum Hash {
    UNSIGNED(u64),
    SIGNED(i64),
}

pub trait HashableChar {
    fn hash_char(&self) -> Hash;
}

macro_rules! impl_hashable_char {
    ($base_type:ty, $kind:tt $(, $t:ty)*) => {
        impl HashableChar for $base_type {
            #[inline]
            fn hash_char(&self) -> Hash
            {
                Hash::$kind(*self $(as $t)*)
            }
        }

        impl HashableChar for &$base_type {
            #[inline]
            fn hash_char(&self) -> Hash
            {
                Hash::$kind(**self $(as $t)*)
            }
        }
    }
}

impl_hashable_char!(char, UNSIGNED, u32, u64);
impl_hashable_char!(i8, SIGNED, i64);
impl_hashable_char!(i16, SIGNED, i64);
impl_hashable_char!(i32, SIGNED, i64);
impl_hashable_char!(i64, SIGNED, i64);
impl_hashable_char!(u8, UNSIGNED, u64);
impl_hashable_char!(u16, UNSIGNED, u64);
impl_hashable_char!(u32, UNSIGNED, u64);
impl_hashable_char!(u64, UNSIGNED, u64);

/// wrapper around iterator which allows iterating over values instead of references
/// todo is there really no way to achieve this without this hack?
#[derive(Clone)]
pub(crate) struct UnrefIterator<'a, T>
where
    T: Copy,
{
    pub seq: core::slice::Iter<'a, T>,
}

impl<'a, T> Iterator for UnrefIterator<'a, T>
where
    T: Copy,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.seq.next().copied()
    }
}
