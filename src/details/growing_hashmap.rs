use crate::{Hash, HashableChar};

#[derive(Default, Clone)]
struct GrowingHashmapMapElem<ValueType> {
    key: u64,
    value: ValueType,
}

/// specialized hashmap to store user provided types
/// this implementation relies on a couple of base assumptions in order to simplify the implementation
/// - the hashmap does not have an upper limit of included items
/// - the default value for the `ValueType` can be used as a dummy value to indicate an empty cell
/// - elements can't be removed
/// - only allocates memory on first write access.
///   This improves performance for hashmaps that are never written to
pub struct GrowingHashmap<ValueType> {
    used: i32,
    fill: i32,
    mask: i32,
    map: Option<Vec<GrowingHashmapMapElem<ValueType>>>,
}

impl<ValueType> Default for GrowingHashmap<ValueType>
where
    ValueType: Default + Clone + Eq,
{
    #[inline]
    fn default() -> Self {
        Self {
            used: 0,
            fill: 0,
            mask: -1,
            map: None,
        }
    }
}

impl<ValueType> GrowingHashmap<ValueType>
where
    ValueType: Default + Clone + Eq + Copy,
{
    #[allow(dead_code)]
    pub const fn size(&self) -> i32 {
        self.used
    }

    #[allow(dead_code)]
    pub const fn capacity(&self) -> i32 {
        self.mask + 1
    }

    #[allow(dead_code)]
    pub const fn empty(&self) -> bool {
        self.used == 0
    }

    pub fn get(&self, key: u64) -> ValueType {
        self.map
            .as_ref()
            .map_or_else(|| Default::default(), |map| map[self.lookup(key)].value)
    }

    pub fn get_mut(&mut self, key: u64) -> &mut ValueType {
        if self.map.is_none() {
            self.allocate();
        }

        let mut i = self.lookup(key);
        if self
            .map
            .as_ref()
            .expect("map should have been created above")[i]
            .value
            == Default::default()
        {
            self.fill += 1;
            // resize when 2/3 full
            if self.fill * 3 >= (self.mask + 1) * 2 {
                self.grow((self.used + 1) * 2);
                i = self.lookup(key);
            }

            self.used += 1;
        }

        let elem = &mut self
            .map
            .as_mut()
            .expect("map should have been created above")[i];
        elem.key = key;
        &mut elem.value
    }

    fn allocate(&mut self) {
        self.mask = 8 - 1;
        self.map = Some(vec![GrowingHashmapMapElem::default(); 8]);
    }

    /// lookup key inside the hashmap using a similar collision resolution
    /// strategy to `CPython` and `Ruby`
    fn lookup(&self, key: u64) -> usize {
        let hash = key;
        let mut i = (hash & self.mask as u64) as usize;

        let map = self
            .map
            .as_ref()
            .expect("callers have to ensure map is allocated");

        if map[i].value == Default::default() || map[i].key == key {
            return i;
        }

        let mut perturb = key;
        loop {
            i = (i * 5 + perturb as usize + 1) % 128;

            if map[i].value == Default::default() || map[i].key == key {
                return i;
            }

            perturb >>= 5;
        }
    }

    fn grow(&mut self, min_used: i32) {
        let mut new_size = self.mask + 1;
        while new_size <= min_used {
            new_size <<= 1;
        }

        let mut new_map = vec![GrowingHashmapMapElem::<ValueType>::default(); new_size as usize];

        self.fill = self.used;
        self.mask = new_size - 1;

        for elem in self
            .map
            .as_ref()
            .expect("callers have to ensure map is allocated")
        {
            if elem.value != Default::default() {
                let j = self.lookup(elem.key);
                let new_elem = &mut new_map[j];
                new_elem.key = elem.key;
                new_elem.value = elem.value;
                self.used -= 1;
                if self.used == 0 {
                    break;
                }
            }
        }

        self.used = self.fill;
        self.map = Some(new_map);
    }
}

pub struct HybridGrowingHashmap<ValueType> {
    // todo in theory we have a fixed keytype here and so we wouldn't need both
    // an unsigned and signed map. In Practice this probably doesn't matter all that much
    pub map_unsigned: GrowingHashmap<ValueType>,
    pub map_signed: GrowingHashmap<ValueType>,
    pub extended_ascii: [ValueType; 256],
}

impl<ValueType> HybridGrowingHashmap<ValueType>
where
    ValueType: Default + Clone + Copy + Eq,
{
    // right now this can't be used since rust fails to elide the memcpy
    // on return
    /*pub fn new() -> Self {
        HybridGrowingHashmap {
            map_unsigned: GrowingHashmap::default(),
            map_signed: GrowingHashmap::default(),
            extended_ascii: [Default::default(); 256],
        }
    }*/

    pub fn get<CharT>(&self, key: CharT) -> ValueType
    where
        CharT: HashableChar,
    {
        match key.hash_char() {
            Hash::SIGNED(value) => {
                if value < 0 {
                    self.map_signed.get(u64::from_ne_bytes(value.to_ne_bytes()))
                } else if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    self.extended_ascii[usize::from(val_u8)]
                } else {
                    self.map_unsigned
                        .get(u64::from_ne_bytes(value.to_ne_bytes()))
                }
            }
            Hash::UNSIGNED(value) => {
                if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    self.extended_ascii[usize::from(val_u8)]
                } else {
                    self.map_unsigned.get(value)
                }
            }
        }
    }

    pub fn get_mut<CharT>(&mut self, key: CharT) -> &mut ValueType
    where
        CharT: HashableChar,
    {
        match key.hash_char() {
            Hash::SIGNED(value) => {
                if value < 0 {
                    self.map_signed
                        .get_mut(u64::from_ne_bytes(value.to_ne_bytes()))
                } else if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    &mut self.extended_ascii[usize::from(val_u8)]
                } else {
                    self.map_unsigned
                        .get_mut(u64::from_ne_bytes(value.to_ne_bytes()))
                }
            }
            Hash::UNSIGNED(value) => {
                if value <= 255 {
                    let val_u8 = u8::try_from(value).expect("we check the bounds above");
                    &mut self.extended_ascii[usize::from(val_u8)]
                } else {
                    self.map_unsigned.get_mut(value)
                }
            }
        }
    }
}
