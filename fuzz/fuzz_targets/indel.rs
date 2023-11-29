#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use rapidfuzz::distance::indel;

#[derive(Arbitrary, Debug)]
pub struct Texts {
    pub s1: String,
    pub s2: String,
}

fn fuzz(texts: Texts) {
    indel::distance(texts.s1.chars(), texts.s2.chars(), None, None).expect("does not return None");

    indel::BatchComparator::new(texts.s1.chars())
        .distance(texts.s2.chars(), None, None)
        .expect("does not return None");
}

fuzz_target!(|texts: Texts| {
    fuzz(texts);
});
