Changelog
---------

[0.3.2] - 2023-11-29
^^^^^^^^^^^^^^^^^^^^
Fixed
~~~~~
* fixed crash inside hashmap grow function leading to a crash in the
  Damerau-Levenshtein implementation

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

