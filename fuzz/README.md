# rapidfuzz-fuzz

This directory contains fuzzers which can be used to automatically identify faults present in Boa. All the fuzzers in
this directory are [grammar-aware](https://www.fuzzingbook.org/html/Grammars.html) (based on
[Arbitrary](https://docs.rs/arbitrary/latest/arbitrary/)) and coverage-guided.

You can run any fuzzer you wish with the following command (replacing `your-fuzzer` with a fuzzer availble in
fuzz_targets, e.g. `damerau-levenshtein`):

```bash
cargo fuzz run -s none your-fuzzer
```

Note that you may wish to use a different sanitizer option (`-s`) according to what kind of issue you're looking for.
Refer to the [cargo-fuzz book](https://rust-fuzz.github.io/book/cargo-fuzz.html) for details on how to select a
sanitizer and other flags.
