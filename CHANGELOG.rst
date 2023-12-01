Changelog
---------

[0.4.0] - 2023-12-01
^^^^^^^^^^^^^^^^^^^^
Changed
~~~~~~~
* rewrite of function signatures to reduce boilerplate

  * return type now automatically deduced, so no more unwrapping needed
    when ``score_cutoff`` is not used
  * optional arguments now in Arg structs uisng the builder pattern to reduce amount
    of extra arguments
  * extra overload ``*_with_args`` for a variant with args, while the default version accepts
    only two sequences

  The signatures is expected to largely stay this way for the foreseeable future.

[0.3.2] - 2023-11-29
^^^^^^^^^^^^^^^^^^^^
Fixed
~~~~~
* fixed crash inside hashmap grow function leading to a crash in the
  Damerau-Levenshtein implementation
* fixed incorrect flagging of similar characters in Jaro similarity
* fixed wraparound in Longest Common Subsequence

[0.3.1] - 2023-11-29
^^^^^^^^^^^^^^^^^^^^
Fixed
~~~~~
* fixed crash inside hashmap lookup function leading to a crash in the
  Damerau-Levenshtein implementation

[0.3.0] - 2023-11-27
^^^^^^^^^^^^^^^^^^^^
Previous versions only existed for testing purposed years ago. This is a complete
rewrite porting a subset of the features provided in the C++ implementation of
rapidfuzz. The remaining features will be added in later releases.

Added
~~~~~
* added implementations of the following string metrics:

  * Levenshtein distance
  * Damerau-Levenshtein distance
  * Hamming distance
  * Longest common subsequence
  * Indel distance
  * Optimal string alignment distance
  * Postfix similarity
  * Prefix similarity
  * Jaro similarity
  * Jaro-Winkler similarity

