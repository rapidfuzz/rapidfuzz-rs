use crate::details::intrinsics::ceil_div_usize;
use crate::details::matrix::BitMatrix;
use crate::{Hash, HashableChar};

#[derive(Clone, Copy, Default)]
struct BitvectorHashmapMapElem {
    key: u64,
    value: u64,
}

/// specialized hashmap to store bitvectors
/// this implementation relies on a couple of base assumptions in order to simplify the implementation
/// - the hashmap includes at most 64 different items
/// - since bitvectors are only in use when at least one bit is set, 0 can be used to indicate an unused element
/// - elements are never explicitly removed. When changing a sliding window over a string, shifting the corresponding
///   bits would eventually be 0 -> removed the element
/// - works with u64 keys. The caller has to ensure these have no collisions when using e.g. a mixture of u64 and i64 elements
///   this can be done e.g. by using two hashmaps one for values < 0 and one for values >= 0
#[derive(Clone)]
pub struct BitvectorHashmap {
    map: [BitvectorHashmapMapElem; 128],
}

impl Default for BitvectorHashmap {
    #[inline]
    fn default() -> Self {
        Self {
            map: [BitvectorHashmapMapElem::default(); 128],
        }
    }
}

impl BitvectorHashmap {
    pub const fn get(&self, key: u64) -> u64 {
        self.map[self.lookup(key)].value
    }

    pub fn get_mut(&mut self, key: u64) -> &mut u64 {
        let i = self.lookup(key);
        let elem = &mut self.map[i];
        elem.key = key;
        &mut elem.value
    }

    /// lookup key inside the hashmap using a similar collision resolution
    /// strategy to `CPython` and `Ruby`
    const fn lookup(&self, key: u64) -> usize {
        let mut i = (key % 128) as usize;

        if self.map[i].value == 0 || self.map[i].key == key {
            return i;
        }

        let mut perturb = key;
        loop {
            i = (i * 5 + perturb as usize + 1) % 128;

            if self.map[i].value == 0 || self.map[i].key == key {
                return i;
            }

            perturb >>= 5;
        }
    }
}

pub struct PatternMatchVector {
    pub extended_ascii: [u64; 256],
    pub map_unsigned: Option<BitvectorHashmap>,
    pub map_signed: Option<BitvectorHashmap>,
}

pub trait BitVectorInterface {
    fn get<CharT>(&self, block: usize, key: CharT) -> u64
    where
        CharT: HashableChar;

    fn size(&self) -> usize;
}

impl Default for PatternMatchVector {
    fn default() -> Self {
        Self {
            map_unsigned: None,
            map_signed: None,
            extended_ascii: [0; 256],
        }
    }
}

impl PatternMatchVector {
    pub fn insert<Iter1, CharT>(&mut self, s1: Iter1)
    where
        Iter1: Iterator<Item = CharT>,
        CharT: HashableChar,
    {
        let mut mask: u64 = 1;
        for ch in s1 {
            self.insert_mask(ch, mask);
            mask <<= 1;
        }
    }

    fn insert_mask<CharT>(&mut self, key: CharT, mask: u64)
    where
        CharT: HashableChar,
    {
        match key.hash_char() {
            Hash::SIGNED(value) => {
                if value < 0 {
                    if self.map_signed.is_none() {
                        self.map_signed = Some(BitvectorHashmap::default());
                    }
                    let item = self
                        .map_signed
                        .as_mut()
                        .expect("map should have been created above")
                        .get_mut(u64::from_ne_bytes(value.to_ne_bytes()));
                    *item |= mask;
                } else if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    let item = &mut self.extended_ascii[usize::from(val_u8)];
                    *item |= mask;
                } else {
                    if self.map_unsigned.is_none() {
                        self.map_unsigned = Some(BitvectorHashmap::default());
                    }
                    let item = self
                        .map_unsigned
                        .as_mut()
                        .expect("map should have been created above")
                        .get_mut(u64::from_ne_bytes(value.to_ne_bytes()));
                    *item |= mask;
                }
            }
            Hash::UNSIGNED(value) => {
                if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    let item = &mut self.extended_ascii[usize::from(val_u8)];
                    *item |= mask;
                } else {
                    if self.map_unsigned.is_none() {
                        self.map_unsigned = Some(BitvectorHashmap::default());
                    }
                    let item = self
                        .map_unsigned
                        .as_mut()
                        .expect("map should have been created above")
                        .get_mut(value);
                    *item |= mask;
                }
            }
        }
    }
}

impl BitVectorInterface for PatternMatchVector {
    fn get<CharT>(&self, block: usize, key: CharT) -> u64
    where
        CharT: HashableChar,
    {
        debug_assert!(block == 0);
        match key.hash_char() {
            Hash::SIGNED(value) => {
                if value < 0 {
                    self.map_signed
                        .as_ref()
                        .map_or(0, |map| map.get(u64::from_ne_bytes(value.to_ne_bytes())))
                } else if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    self.extended_ascii[usize::from(val_u8)]
                } else {
                    self.map_unsigned
                        .as_ref()
                        .map_or(0, |map| map.get(u64::from_ne_bytes(value.to_ne_bytes())))
                }
            }
            Hash::UNSIGNED(value) => {
                if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    self.extended_ascii[usize::from(val_u8)]
                } else {
                    self.map_unsigned.as_ref().map_or(0, |map| map.get(value))
                }
            }
        }
    }

    fn size(&self) -> usize {
        1
    }
}

#[derive(Clone)]
pub struct BlockPatternMatchVector {
    pub block_count: usize,
    pub map_unsigned: Option<Vec<BitvectorHashmap>>,
    pub map_signed: Option<Vec<BitvectorHashmap>>,
    pub extended_ascii: BitMatrix<u64>,
}

impl BlockPatternMatchVector {
    pub fn new(str_len: usize) -> Self {
        let block_count = ceil_div_usize(str_len, 64);
        Self {
            block_count,
            map_unsigned: None,
            map_signed: None,
            extended_ascii: BitMatrix::<u64>::new(256, block_count, 0),
        }
    }

    pub fn insert<Iter1, CharT>(&mut self, s1: Iter1)
    where
        Iter1: Iterator<Item = CharT>,
        CharT: HashableChar,
    {
        let mut mask: u64 = 1;
        for (i, ch) in s1.enumerate() {
            let block = i / 64;
            self.insert_mask(block, ch, mask);
            mask = mask.rotate_left(1);
        }
    }

    fn insert_mask<CharT>(&mut self, block: usize, key: CharT, mask: u64)
    where
        CharT: HashableChar,
    {
        debug_assert!(block < self.size());

        match key.hash_char() {
            Hash::SIGNED(value) => {
                if value < 0 {
                    if self.map_signed.is_none() {
                        self.map_signed = Some(vec![BitvectorHashmap::default(); self.block_count]);
                    }
                    let item = self
                        .map_signed
                        .as_mut()
                        .expect("map should have been created above")[block]
                        .get_mut(u64::from_ne_bytes(value.to_ne_bytes()));
                    *item |= mask;
                } else if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    let item = self.extended_ascii.get_mut(val_u8.into(), block);
                    *item |= mask;
                } else {
                    if self.map_unsigned.is_none() {
                        self.map_unsigned =
                            Some(vec![BitvectorHashmap::default(); self.block_count]);
                    }
                    let item = self
                        .map_unsigned
                        .as_mut()
                        .expect("map should have been created above")[block]
                        .get_mut(u64::from_ne_bytes(value.to_ne_bytes()));
                    *item |= mask;
                }
            }
            Hash::UNSIGNED(value) => {
                if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    let item = self.extended_ascii.get_mut(val_u8.into(), block);
                    *item |= mask;
                } else {
                    if self.map_unsigned.is_none() {
                        self.map_unsigned =
                            Some(vec![BitvectorHashmap::default(); self.block_count]);
                    }
                    let item = self
                        .map_unsigned
                        .as_mut()
                        .expect("map should have been created above")[block]
                        .get_mut(value);
                    *item |= mask;
                }
            }
        }
    }
}

impl BitVectorInterface for BlockPatternMatchVector {
    fn get<CharT>(&self, block: usize, key: CharT) -> u64
    where
        CharT: HashableChar,
    {
        debug_assert!(block < self.size());

        match key.hash_char() {
            Hash::SIGNED(value) => {
                if value < 0 {
                    self.map_signed.as_ref().map_or(0, |map| {
                        map[block].get(u64::from_ne_bytes(value.to_ne_bytes()))
                    })
                } else if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    *self.extended_ascii.get(val_u8.into(), block)
                } else {
                    self.map_unsigned.as_ref().map_or(0, |map| {
                        map[block].get(u64::from_ne_bytes(value.to_ne_bytes()))
                    })
                }
            }
            Hash::UNSIGNED(value) => {
                if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    *self.extended_ascii.get(val_u8.into(), block)
                } else {
                    self.map_unsigned
                        .as_ref()
                        .map_or(0, |map| map[block].get(value))
                }
            }
        }
    }

    fn size(&self) -> usize {
        self.block_count
    }
}
