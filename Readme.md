<h1 align="center">
<img src="https://raw.githubusercontent.com/maxbachmann/rapidfuzz/master/docs/img/RapidFuzz.svg?sanitize=true" alt="RapidFuzz" width="400">
</h1>
<h4 align="center">Rapid fuzzy string matching in Rust using the Levenshtein Distance</h4>

<p align="center">
  <a href="https://github.com/maxbachmann/rapidfuzz-rs/actions">
    <img src="https://github.com/maxbachmann/rapidfuzz-rs/workflows/Rust/badge.svg"
         alt="Continous Integration">
  </a>
  <a href="https://gitter.im/rapidfuzz/community">
    <img src="https://badges.gitter.im/rapidfuzz/community.svg"
         alt="Gitter chat">
  </a>
  <a href="https://maxbachmann.github.io/rapidfuzz">
    <img src="https://img.shields.io/badge/-documentation-blue"
         alt="Documentation">
  </a>
  <a href="https://img.shields.io/crates/l/rapidfuzz.svg">
    <img src="https://img.shields.io/crates/l/rapidfuzz.svg"
         alt="license">
  </a>
</p>

<p align="center">
  <a href="#description">Description</a> •
  <a href="#installation">Installation</a> •
  <a href="#usage">Usage</a> •
  <a href="#license">License</a>
</p>

---
## Description

RapidFuzz is a general purpose string matching library with implementations
for Rust, C++ and Python.

### Key Features

- **Diverse String Metrics**: Offers a variety of string metrics
  to suit different use cases. These range from the Levenshtein
  distance for edit-based comparisons to the Jaro-Winkler similarity for
  more nuanced similarity assessments.
- **Optimized for Speed**: The library is designed with performance in mind.
  Each implementation is carefully designed to ensure optimal performance,
  making it suitable for the analysis of large datasets.
- **Easy to use**: The API is designed to be simple to use, while still giving
  the implementation room for optimization.

## Installation

The installation is as simple as:
```console
$ cargo add rapidfuzz
```

## Usage

The following examples show the usage with the Levenshtein distance. Other metrics
can be found in the [fuzz](https://docs.rs/rapidfuzz/latest/rapidfuzz/fuzz/index.html) and [distance](https://docs.rs/rapidfuzz/latest/rapidfuzz/distance/index.html) modules.

```rust
use rapidfuzz::distance::levenshtein;

// Perform a simple comparision using he levenshtein distance
assert_eq!(
    Some(3),
    levenshtein::distance("kitten".chars(), "sitting".chars(), None, None, None)
);

// If you are sure the input strings are ASCII only it's usually faster to operate on bytes
assert_eq!(
    Some(3),
    levenshtein::distance("kitten".bytes(), "sitting".bytes(), None, None, None)
);

// You can provide a score_cutoff value to filter out strings with distance that is worse than
// the score_cutoff
assert_eq!(
    None,
    levenshtein::distance("kitten".chars(), "sitting".chars(), None, 2, None)
);

// You can provide a score_hint to tell the implementation about the expected score.
// This can be used to select a more performant implementation internally, but might cause
// a slowdown in cases where the distance is actually worse than the score_hint
assert_eq!(
    Some(3),
    levenshtein::distance("kitten".chars(), "sitting".chars(), None, None, Some(3))
);

// When comparing a single string to multiple strings you can use the provided `BatchComparators`.
// These can cache part of the calculation which can provide significant speedups
let scorer = levenshtein::BatchComparator::new("kitten".chars(), None);
assert_eq!(Some(3), scorer.distance("sitting".chars(), None, None));
assert_eq!(Some(0), scorer.distance("kitten".chars(), None, None));
```


## License
Licensed under either of [Apache License, Version
2.0](https://github.com/maxbachmann/rapidfuzz-rs/blob/main/LICENSE-APACHE) or [MIT License](https://github.com/maxbachmann/rapidfuzz-rs/blob/main/LICENSE-MIT) at your option.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in RapidFuzz by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

