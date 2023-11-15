use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::{distributions::Alphanumeric, Rng};

use rapidfuzz::distance::{
    damerau_levenshtein_distance, hamming_distance, levenshtein_distance, CachedDamerauLevenshtein,
    CachedLevenshtein, LevenshteinWeightTable,
};

use std::str::Bytes;
use strsim::{generic_damerau_levenshtein, generic_levenshtein};

fn generate(len: usize) -> String {
    rand::thread_rng()
        .sample_iter(&Alphanumeric)
        .take(len)
        .map(char::from)
        .collect()
}

struct StringWrapper<'a>(&'a str);

impl<'a, 'b> IntoIterator for &'a StringWrapper<'b> {
    type Item = u8;
    type IntoIter = Bytes<'b>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.bytes()
    }
}

fn benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Levenshtein");

    for i in [4, 6, 8, 10, 12, 16, 32, 64, 128, 256].iter() {
        let s1 = generate(*i);
        let s2 = generate(*i);

        group.bench_with_input(BenchmarkId::new("rapidfuzz", i), &(&s1, &s2), |b, val| {
            b.iter(|| {
                black_box(levenshtein_distance(
                    val.0.bytes(),
                    val.1.bytes(),
                    None,
                    None,
                    None,
                ));
            })
        });
        group.bench_with_input(BenchmarkId::new("strsim", i), &(&s1, &s2), |b, val| {
            b.iter(|| {
                black_box(generic_levenshtein(
                    &StringWrapper(val.0),
                    &StringWrapper(val.1),
                ));
            })
        });

        let cached = CachedLevenshtein::new(s1.bytes(), None);
        group.bench_with_input(
            BenchmarkId::new("cached_rapidfuzz", i),
            &(&cached, &s2),
            |b, val| {
                b.iter(|| {
                    black_box(cached.distance(val.1.bytes(), None, None));
                })
            },
        );
    }

    group.finish();

    group = c.benchmark_group("DamerauLevenshtein");

    for i in [4, 6, 8, 10, 12, 16, 32, 64, 128, 256].iter() {
        let s1 = generate(*i);
        let s2 = generate(*i);

        group.bench_with_input(BenchmarkId::new("rapidfuzz", i), &(&s1, &s2), |b, val| {
            b.iter(|| {
                black_box(damerau_levenshtein_distance(
                    val.0.bytes(),
                    val.1.bytes(),
                    None,
                    None,
                ));
            })
        });
        let (x, y): (Vec<_>, Vec<_>) = (s1.bytes().collect(), s2.bytes().collect());
        group.bench_with_input(BenchmarkId::new("strsim", i), &(&x, &y), |b, val| {
            b.iter(|| {
                black_box(generic_damerau_levenshtein(
                    val.0.as_slice(),
                    val.1.as_slice(),
                ));
            })
        });
    }

    group.finish();
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
