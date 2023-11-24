use crate::details::common::{norm_sim_to_norm_dist, HashableChar};

pub trait DistanceMetricUsize {
    fn maximum(&self, len1: usize, len2: usize) -> usize;

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        mut score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);
        if score_cutoff > maximum {
            return 0;
        }

        score_hint = score_hint.min(score_cutoff);
        let cutoff_distance = maximum - score_cutoff;
        let hint_distance = maximum - score_hint;
        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance);
        let sim = maximum - dist;
        if sim >= score_cutoff {
            sim
        } else {
            0
        }
    }
}

pub trait SimilarityMetricUsize {
    fn maximum(&self, len1: usize, len2: usize) -> usize;

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_similarity = if maximum >= score_cutoff {
            maximum - score_cutoff
        } else {
            0
        };
        let hint_similarity = if maximum >= score_hint {
            maximum - score_hint
        } else {
            0
        };

        let sim = self._similarity(s1, len1, s2, len2, cutoff_similarity, hint_similarity);
        let dist = maximum - sim;

        if dist <= score_cutoff {
            dist
        } else {
            score_cutoff + 1
        }
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: usize,
        score_hint: usize,
    ) -> usize
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
}

pub trait NormalizedMetricUsize {
    fn _normalized_distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
}

impl<T: DistanceMetricUsize> NormalizedMetricUsize for T {
    fn _normalized_distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_distance = (maximum as f64 * score_cutoff).ceil() as usize;
        let hint_distance = (maximum as f64 * score_hint).ceil() as usize;

        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance);
        let norm_dist = if maximum == 0 {
            0.0
        } else {
            dist as f64 / maximum as f64
        };
        if norm_dist <= score_cutoff {
            norm_dist
        } else {
            1.0
        }
    }

    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let cutoff_score = norm_sim_to_norm_dist(score_cutoff);
        let hint_score = norm_sim_to_norm_dist(score_hint);

        let norm_dist = self._normalized_distance(s1, len1, s2, len2, cutoff_score, hint_score);
        let norm_sim = 1.0 - norm_dist;

        if norm_sim >= score_cutoff {
            norm_sim
        } else {
            0.0
        }
    }
}

// todo how to deduplicate this? Right now NormalizedMetricUsize2
// could be placed inside SimilarityMetricUsize since it's duplicated anyways
pub trait NormalizedMetricUsize2 {
    fn _normalized_distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
}

impl<T: SimilarityMetricUsize> NormalizedMetricUsize2 for T {
    fn _normalized_distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_distance = (maximum as f64 * score_cutoff).ceil() as usize;
        let hint_distance = (maximum as f64 * score_hint).ceil() as usize;

        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance);
        let norm_dist = if maximum == 0 {
            0.0
        } else {
            dist as f64 / maximum as f64
        };
        if norm_dist <= score_cutoff {
            norm_dist
        } else {
            1.0
        }
    }

    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> f64
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let cutoff_score = norm_sim_to_norm_dist(score_cutoff);
        let hint_score = norm_sim_to_norm_dist(score_hint);

        let norm_dist = self._normalized_distance(s1, len1, s2, len2, cutoff_score, hint_score);
        let norm_sim = 1.0 - norm_dist;

        if norm_sim >= score_cutoff {
            norm_sim
        } else {
            0.0
        }
    }
}

pub trait SimilarityMetricf64 {
    fn maximum(&self, len1: usize, len2: usize) -> f64;

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_similarity = if maximum >= score_cutoff {
            maximum - score_cutoff
        } else {
            0.0
        };
        let hint_similarity = if maximum >= score_hint {
            maximum - score_hint
        } else {
            0.0
        };

        let sim = self._similarity(s1, len1, s2, len2, cutoff_similarity, hint_similarity)?;
        let dist = maximum - sim;

        if dist <= score_cutoff {
            Some(dist)
        } else {
            None
        }
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
}

pub trait NormalizedMetricf64 {
    fn _normalized_distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy;
}

impl<T: SimilarityMetricf64> NormalizedMetricf64 for T {
    fn _normalized_distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_distance = maximum * score_cutoff;
        let hint_distance = maximum * score_hint;

        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance)?;
        let norm_dist = if maximum > 0.0 { dist / maximum } else { 0.0 };
        if norm_dist <= score_cutoff {
            Some(norm_dist)
        } else {
            None
        }
    }

    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: f64,
        score_hint: f64,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let cutoff_score = norm_sim_to_norm_dist(score_cutoff);
        let hint_score = norm_sim_to_norm_dist(score_hint);

        let norm_dist = self._normalized_distance(s1, len1, s2, len2, cutoff_score, hint_score)?;
        let norm_sim = 1.0 - norm_dist;

        if norm_sim >= score_cutoff {
            Some(norm_sim)
        } else {
            None
        }
    }
}
