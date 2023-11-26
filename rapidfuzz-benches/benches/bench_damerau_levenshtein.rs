use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::{distributions::Alphanumeric, Rng};

use rapidfuzz::distance;

use std::str::Bytes;

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
    let mut group = c.benchmark_group("Damerau Levenshtein");

    for i in (2..128).step_by(2) {
        let s1 = generate(i);
        let s2 = generate(i);

        group.bench_with_input(BenchmarkId::new("rapidfuzz", i), &(&s1, &s2), |b, val| {
            b.iter(|| {
                black_box(distance::damerau_levenshtein::distance(
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
                black_box(strsim::generic_damerau_levenshtein(
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
