use crate::details::common::norm_sim_to_norm_dist;
use crate::HashableChar;

pub trait MetricUsize {
    fn maximum(&self, len1: usize, len2: usize) -> usize;

    // todo rename if we convert everything to use the builder pattern
    fn distance_<Iter1, Iter2>(
        &self,
        s1: Iter1,
        s2: Iter2,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let s1_iter = s1.into_iter();
        let s2_iter = s2.into_iter();
        self._distance(
            s1_iter.clone(),
            s1_iter.count(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff,
            score_hint,
        )
    }

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
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

    // todo rename if we convert everything to use the builder pattern
    fn similarity_<Iter1, Iter2>(
        &self,
        s1: Iter1,
        s2: Iter2,
        score_cutoff: Option<usize>,
        score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let s1_iter = s1.into_iter();
        let s2_iter = s2.into_iter();
        self._similarity(
            s1_iter.clone(),
            s1_iter.count(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff,
            score_hint,
        )
    }

    fn _similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<usize>,
        mut score_hint: Option<usize>,
    ) -> Option<usize>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);
        if let Some(cutoff) = score_cutoff {
            if cutoff > maximum {
                return None;
            }

            if let Some(hint) = score_hint {
                score_hint = Some(hint.min(cutoff));
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
        Some(sim)
    }

    // todo rename if we convert everything to use the builder pattern
    fn normalized_distance_<Iter1, Iter2>(
        &self,
        s1: Iter1,
        s2: Iter2,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let s1_iter = s1.into_iter();
        let s2_iter = s2.into_iter();
        self._normalized_distance(
            s1_iter.clone(),
            s1_iter.count(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff,
            score_hint,
        )
    }

    fn _normalized_distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        mut score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);

        let cutoff_distance;
        if let Some(mut cutoff) = score_cutoff {
            cutoff = cutoff.clamp(0.0, 1.0);
            score_cutoff = Some(cutoff);
            cutoff_distance = Some((maximum as f64 * cutoff).ceil() as usize);
        } else {
            cutoff_distance = None;
        }

        let hint_distance;
        if let Some(mut cutoff) = score_hint {
            cutoff = cutoff.clamp(0.0, 1.0);
            hint_distance = Some((maximum as f64 * cutoff).ceil() as usize);
        } else {
            hint_distance = None;
        }

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

    // todo rename if we convert everything to use the builder pattern
    fn normalized_similarity_<Iter1, Iter2>(
        &self,
        s1: Iter1,
        s2: Iter2,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: IntoIterator,
        Iter1::IntoIter: DoubleEndedIterator + Clone,
        Iter2: IntoIterator,
        Iter2::IntoIter: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let s1_iter = s1.into_iter();
        let s2_iter = s2.into_iter();
        self._normalized_similarity(
            s1_iter.clone(),
            s1_iter.count(),
            s2_iter.clone(),
            s2_iter.count(),
            score_cutoff,
            score_hint,
        )
    }

    fn _normalized_similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
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

pub trait Metricf64 {
    fn maximum(&self, len1: usize, len2: usize) -> f64;

    fn _distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
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

    fn _similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        mut score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
    {
        let maximum = self.maximum(len1, len2);
        if let Some(cutoff) = score_cutoff {
            if cutoff > maximum {
                return None;
            }

            if let Some(hint) = score_hint {
                score_hint = Some(hint.min(cutoff));
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
        Some(sim)
    }

    fn _normalized_distance<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
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

    fn _normalized_similarity<Iter1, Iter2>(
        &self,
        s1: Iter1,
        len1: usize,
        s2: Iter2,
        len2: usize,
        score_cutoff: Option<f64>,
        score_hint: Option<f64>,
    ) -> Option<f64>
    where
        Iter1: DoubleEndedIterator + Clone,
        Iter2: DoubleEndedIterator + Clone,
        Iter1::Item: PartialEq<Iter2::Item> + HashableChar + Copy,
        Iter2::Item: PartialEq<Iter1::Item> + HashableChar + Copy,
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
