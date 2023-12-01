#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use rapidfuzz::distance::jaro_winkler;

#[derive(Arbitrary, Debug)]
pub struct Texts {
    pub s1: String,
    pub s2: String,
}

fn fuzz(texts: Texts) {
    jaro_winkler::distance(texts.s1.chars(), texts.s2.chars());

    jaro_winkler::BatchComparator::new(texts.s1.chars()).distance(texts.s2.chars());
}

fuzz_target!(|texts: Texts| {
    fuzz(texts);
});
