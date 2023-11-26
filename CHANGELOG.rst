Changelog
---------

[0.3.0] - 2023-
^^^^^^^^^^^^^^^^^^^^^
Previous versions only existed for testing purposed years ago. This is a complete
rewrite porting a subset of the features provided in the C++ implementation of
rapidfuzz. The remaining features will be added in later releases.

Added
~~~~~~~
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
