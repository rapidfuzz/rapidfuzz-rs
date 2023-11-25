use crate::details::common::{norm_sim_to_norm_dist, HashableChar};

pub trait DistanceMetricUsize {
    fn maximum(&self, len1: usize, len2: usize) -> usize;

    fn _distance<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
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
        score_cutoff: Option<usize>,
        mut score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);
        if let Some(cutoff) = score_cutoff {
            if cutoff > maximum {
                return None;
            }

            if let Some(hint) = score_hint {
                score_hint = Some(hint.min(cutoff))
            }
        }

        let cutoff_distance = score_cutoff.map(|x| maximum - x);
        let hint_distance = score_hint.map(|x| maximum - x);
        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance)?;
        let sim = maximum - dist;
        if let Some(cutoff) = score_cutoff {
            if sim < cutoff {
                return None;
            }
        }
        Some(dist)
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
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_similarity = score_cutoff.map(|x| if maximum >= x { maximum - x } else { 0 });
        let hint_similarity = score_hint.map(|x| if maximum >= x { maximum - x } else { 0 });

        let sim = self._similarity(s1, len1, s2, len2, cutoff_similarity, hint_similarity)?;
        let dist = maximum - sim;

        if let Some(cutoff) = score_cutoff {
            if dist > cutoff {
                return None;
            }
        }
        Some(dist)
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_distance = score_cutoff.map(|x| (maximum as f64 * x).ceil() as usize);
        let hint_distance = score_hint.map(|x| (maximum as f64 * x).ceil() as usize);

        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance)?;
        let norm_dist = if maximum == 0 {
            0.0
        } else {
            dist as f64 / maximum as f64
        };
        if let Some(cutoff) = score_cutoff {
            if norm_dist > cutoff {
                return None;
            }
        }
        Some(norm_dist)
    }

    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let cutoff_score = score_cutoff.map(norm_sim_to_norm_dist);
        let hint_score = score_hint.map(norm_sim_to_norm_dist);

        let norm_dist = self._normalized_distance(s1, len1, s2, len2, cutoff_score, hint_score)?;
        let norm_sim = 1.0 - norm_dist;

        if let Some(cutoff) = score_cutoff {
            if norm_sim < cutoff {
                return None;
            }
        }
        Some(norm_sim)
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_distance = score_cutoff.map(|x| (maximum as f64 * x).ceil() as usize);
        let hint_distance = score_hint.map(|x| (maximum as f64 * x).ceil() as usize);

        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance)?;
        let norm_dist = if maximum == 0 {
            0.0
        } else {
            dist as f64 / maximum as f64
        };
        if let Some(cutoff) = score_cutoff {
            if norm_dist > cutoff {
                return None;
            }
        }
        Some(norm_dist)
    }

    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let cutoff_score = score_cutoff.map(norm_sim_to_norm_dist);
        let hint_score = score_hint.map(norm_sim_to_norm_dist);

        let norm_dist = self._normalized_distance(s1, len1, s2, len2, cutoff_score, hint_score)?;
        let norm_sim = 1.0 - norm_dist;

        if let Some(cutoff) = score_cutoff {
            if norm_sim < cutoff {
                return None;
            }
        }
        Some(norm_sim)
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_similarity = score_cutoff.map(|x| if maximum >= x { maximum - x } else { 0.0 });
        let hint_similarity = score_hint.map(|x| if maximum >= x { maximum - x } else { 0.0 });

        let sim = self._similarity(s1, len1, s2, len2, cutoff_similarity, hint_similarity)?;
        let dist = maximum - sim;

        if let Some(cutoff) = score_cutoff {
            if dist > cutoff {
                return None;
            }
        }
        Some(dist)
    }

    fn _similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
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
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_distance = score_cutoff.map(|x| maximum * x);
        let hint_distance = score_hint.map(|x| maximum * x);

        let dist = self._distance(s1, len1, s2, len2, cutoff_distance, hint_distance)?;
        let norm_dist = if maximum > 0.0 { dist / maximum } else { 0.0 };
        if let Some(cutoff) = score_cutoff {
            if norm_dist > cutoff {
                return None;
            }
        }
        Some(norm_dist)
    }

    fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: Iterator<Item = Elem1> + DoubleEndedIterator + Clone,
        Iter2: Iterator<Item = Elem2> + DoubleEndedIterator + Clone,
        Elem1: PartialEq<Elem2> + HashableChar + Copy,
        Elem2: PartialEq<Elem1> + HashableChar + Copy,
    {
        let cutoff_score = score_cutoff.map(norm_sim_to_norm_dist);
        let hint_score = score_hint.map(norm_sim_to_norm_dist);

        let norm_dist = self._normalized_distance(s1, len1, s2, len2, cutoff_score, hint_score)?;
        let norm_sim = 1.0 - norm_dist;

        if let Some(cutoff) = score_cutoff {
            if norm_sim < cutoff {
                return None;
            }
        }
        Some(norm_sim)
    }
}
