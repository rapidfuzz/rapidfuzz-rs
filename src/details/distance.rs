macro_rules! less_than_score_cutoff_similarity {
    ($score_cutoff:expr, f32) => {
        1.0
    };
    ($score_cutoff:expr, f64) => {
        1.0
    };
    ($score_cutoff:expr, $tp:ty) => {
        $score_cutoff + 1
    };
}

// todo maybe some of these could be traits instead?
macro_rules! build_normalized_metric_funcs
{
    ($impl_type:tt, $res_type:ty, $worst_similarity:expr, $worst_distance:expr $(, $v:ident: $t:ty)*) => {
        #[allow(dead_code)]
        pub fn normalized_distance<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            s2: Iter2,
            $($v: $t,)*
            score_cutoff: Option<f64>,
            score_hint: Option<f64>
        ) -> f64
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let len1 = s1_iter.clone().count();
            let len2 = s2_iter.clone().count();
            $impl_type::_normalized_distance(
                s1_iter,
                len1,
                s2_iter,
                len2,
                $($v,)*
                score_cutoff.unwrap_or(1.0),
                score_hint.unwrap_or(1.0)
            )
        }

        pub(crate) fn _normalized_distance<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            len1: usize,
            s2: Iter2,
            len2: usize,
            $($v: $t,)*
            score_cutoff: f64,
            score_hint: f64
        ) -> f64
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let maximum = $impl_type::maximum(s1_iter.clone(), len1, s2_iter.clone(), len2, $($v,)*);

            let cutoff_distance = (maximum as f64 * score_cutoff).ceil() as $res_type;
            let hint_distance = (maximum as f64 * score_hint).ceil() as $res_type;

            let dist = $impl_type::_distance(
                s1_iter,
                len1,
                s2_iter,
                len2,
                $($v,)*
                cutoff_distance,
                hint_distance
            );
            let norm_dist = if maximum != 0 as $res_type {
                dist as f64 / maximum as f64
            } else {
                0.0
            };
            if norm_dist <= score_cutoff {
                norm_dist
            } else {
                1.0
            }
        }

        #[allow(dead_code)]
        pub fn normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            s2: Iter2,
            $($v: $t,)*
            score_cutoff: Option<f64>,
            score_hint: Option<f64>
        ) -> f64
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let len1 = s1_iter.clone().count();
            let len2 = s2_iter.clone().count();
            $impl_type::_normalized_similarity(
                s1_iter,
                len1,
                s2_iter,
                len2,
                $($v,)*
                score_cutoff.unwrap_or(0.0),
                score_hint.unwrap_or(0.0)
            )
        }

        pub(crate) fn _normalized_similarity<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            len1: usize,
            s2: Iter2,
            len2: usize,
            $($v: $t,)*
            score_cutoff: f64,
            score_hint: f64
        ) -> f64
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let cutoff_score = norm_sim_to_norm_dist(score_cutoff);
            let hint_score = norm_sim_to_norm_dist(score_hint);

            let norm_dist = $impl_type::_normalized_distance(
                s1,
                len1,
                s2,
                len2,
                $($v,)*
                cutoff_score,
                hint_score
            );
            let norm_sim = 1.0 - norm_dist;

            if norm_sim >= score_cutoff {
                norm_sim
            } else {
                0.0
            }
        }
    };
}

macro_rules! build_distance_metric_funcs
{
    ($impl_type:tt, $res_type:ty, $worst_similarity:expr, $worst_distance:expr $(, $v:ident: $t:ty)*) => {
        build_normalized_metric_funcs!($impl_type, $res_type, $worst_similarity, $worst_distance $(, $v: $t)*);

        #[allow(dead_code)]
        pub fn distance<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            s2: Iter2,
            $($v: $t,)*
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>
        ) -> $res_type
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let len1 = s1_iter.clone().count();
            let len2 = s2_iter.clone().count();
            $impl_type::_distance(
                s1_iter,
                len1,
                s2_iter,
                len2,
                $($v,)*
                score_cutoff.unwrap_or($worst_distance),
                score_hint.unwrap_or($worst_distance)
            )
        }

        #[allow(dead_code)]
        pub fn similarity<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            s2: Iter2,
            $($v: $t,)*
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>
        ) -> $res_type
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let len1 = s1_iter.clone().count();
            let len2 = s2_iter.clone().count();
            $impl_type::_similarity(
                s1_iter,
                len1,
                s2_iter,
                len2,
                $($v,)*
                score_cutoff.unwrap_or($worst_similarity),
                score_hint.unwrap_or($worst_similarity)
            )
        }

        pub(crate) fn _similarity<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            len1: usize,
            s2: Iter2,
            len2: usize,
            $($v: $t,)*
            score_cutoff: $res_type,
            mut score_hint: $res_type
        ) -> $res_type
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let maximum = $impl_type::maximum(s1_iter.clone(), len1, s2_iter.clone(), len2, $($v,)*);
            if score_cutoff > maximum {
                return 0 as $res_type;
            }

            score_hint = score_hint.min(score_cutoff);
            let cutoff_distance = maximum - score_cutoff;
            let hint_distance = maximum - score_hint;
            let dist = $impl_type::_distance(s1_iter, len1, s2_iter, len2, $($v,)* cutoff_distance, hint_distance);
            let sim = maximum - dist;
            if sim >= score_cutoff {
                sim
            } else {
                0 as $res_type
            }
        }
    };
}

macro_rules! build_similarity_metric_funcs
{
    ($impl_type:tt, $res_type:tt, $worst_similarity:expr, $worst_distance:expr $(, $v:ident: $t:ty)*) => {
        build_normalized_metric_funcs!($impl_type, $res_type, $worst_similarity, $worst_distance $(, $v: $t)*);

        #[allow(dead_code)]
        pub fn distance<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            s2: Iter2,
            $($v: $t,)*
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>
        ) -> $res_type
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let len1 = s1_iter.clone().count();
            let len2 = s2_iter.clone().count();
            $impl_type::_distance(
                s1_iter,
                len1,
                s2_iter,
                len2,
                $($v,)*
                score_cutoff.unwrap_or($worst_distance),
                score_hint.unwrap_or($worst_distance)
            )
        }

        #[allow(dead_code)]
        pub fn similarity<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            s2: Iter2,
            $($v: $t,)*
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>
        ) -> $res_type
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let len1 = s1_iter.clone().count();
            let len2 = s2_iter.clone().count();
            $impl_type::_similarity(
                s1_iter,
                len1,
                s2_iter,
                len2,
                $($v,)*
                score_cutoff.unwrap_or($worst_similarity),
                score_hint.unwrap_or($worst_similarity)
            )
        }

        pub(crate) fn _distance<Iter1, Iter2, Elem1, Elem2>(
            s1: Iter1,
            len1: usize,
            s2: Iter2,
            len2: usize,
            $($v: $t,)*
            score_cutoff: $res_type,
            score_hint: $res_type,
        ) -> $res_type
        where
            Iter1: IntoIterator<Item = Elem1>,
            Iter1::IntoIter: Clone,
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter1 as IntoIterator>::IntoIter: DoubleEndedIterator,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s1_iter = s1.into_iter();
            let s2_iter = s2.into_iter();
            let maximum = $impl_type::maximum(s1_iter.clone(), len1, s2_iter.clone(), len2, $($v,)*);

            let cutoff_similarity = if maximum >= score_cutoff {
                maximum - score_cutoff
            } else {
                $worst_similarity as $res_type
            };
            let hint_similarity = if maximum >= score_hint {
                maximum - score_hint
            } else {
                $worst_similarity as $res_type
            };

            let sim = $impl_type::_similarity(s1_iter, len1, s2_iter, len2, $($v,)* cutoff_similarity, hint_similarity);
            let dist = maximum - sim;

            if dist <= score_cutoff {
                dist
            } else {
                less_than_score_cutoff_similarity!(score_cutoff, $res_type)
            }
        }
    };
}

pub(crate) use build_distance_metric_funcs;
pub(crate) use build_normalized_metric_funcs;
pub(crate) use build_similarity_metric_funcs;
pub(crate) use less_than_score_cutoff_similarity;

macro_rules! build_cached_normalized_metric_funcs {
    ($impl_type:tt, $res_type:ty, $worst_similarity:expr, $worst_distance:expr) => {
        #[allow(dead_code)]
        pub fn normalized_distance<Iter2, Elem2>(
            &self,
            s2: Iter2,
            score_cutoff: Option<f64>,
            score_hint: Option<f64>,
        ) -> f64
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let len2 = s2_iter.clone().count();
            self._normalized_distance(
                s2_iter,
                len2,
                score_cutoff.unwrap_or(1.0),
                score_hint.unwrap_or(1.0),
            )
        }

        pub(crate) fn _normalized_distance<Iter2, Elem2>(
            &self,
            s2: Iter2,
            len2: usize,
            score_cutoff: f64,
            score_hint: f64,
        ) -> f64
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let maximum = self.maximum(s2_iter.clone(), len2);

            let cutoff_distance = (maximum as f64 * score_cutoff).ceil() as $res_type;
            let hint_distance = (maximum as f64 * score_hint).ceil() as $res_type;

            let dist = self._distance(s2_iter, len2, cutoff_distance, hint_distance);
            let norm_dist = if maximum != 0 as $res_type {
                dist as f64 / maximum as f64
            } else {
                0.0
            };
            if norm_dist <= score_cutoff {
                norm_dist
            } else {
                1.0
            }
        }

        #[allow(dead_code)]
        pub fn normalized_similarity<Iter2, Elem2>(
            &self,
            s2: Iter2,
            score_cutoff: Option<f64>,
            score_hint: Option<f64>,
        ) -> f64
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let len2 = s2_iter.clone().count();
            self._normalized_similarity(
                s2_iter,
                len2,
                score_cutoff.unwrap_or(0.0),
                score_hint.unwrap_or(0.0),
            )
        }

        pub(crate) fn _normalized_similarity<Iter2, Elem2>(
            &self,
            s2: Iter2,
            len2: usize,
            score_cutoff: f64,
            score_hint: f64,
        ) -> f64
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let cutoff_score = norm_sim_to_norm_dist(score_cutoff);
            let hint_score = norm_sim_to_norm_dist(score_hint);

            let norm_dist = self._normalized_distance(s2, len2, cutoff_score, hint_score);
            let norm_sim = 1.0 - norm_dist;

            if norm_sim >= score_cutoff {
                norm_sim
            } else {
                0.0
            }
        }
    };
}

macro_rules! build_cached_distance_metric_funcs {
    ($impl_type:tt, $res_type:ty, $worst_similarity:expr, $worst_distance:expr) => {
        build_cached_normalized_metric_funcs!(
            $impl_type,
            $res_type,
            $worst_similarity,
            $worst_distance
        );

        #[allow(dead_code)]
        pub fn distance<Iter2, Elem2>(
            &self,
            s2: Iter2,
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>,
        ) -> $res_type
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let len2 = s2_iter.clone().count();
            self._distance(
                s2_iter,
                len2,
                score_cutoff.unwrap_or($worst_distance),
                score_hint.unwrap_or($worst_distance),
            )
        }

        #[allow(dead_code)]
        pub fn similarity<Iter2, Elem2>(
            &self,
            s2: Iter2,
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>,
        ) -> $res_type
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let len2 = s2_iter.clone().count();
            self._similarity(
                s2_iter,
                len2,
                score_cutoff.unwrap_or($worst_similarity),
                score_hint.unwrap_or($worst_similarity),
            )
        }

        pub(crate) fn _similarity<Iter2, Elem2>(
            &self,
            s2: Iter2,
            len2: usize,
            score_cutoff: $res_type,
            mut score_hint: $res_type,
        ) -> $res_type
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let maximum = self.maximum(s2_iter.clone(), len2);
            if score_cutoff > maximum {
                return 0 as $res_type;
            }

            score_hint = score_hint.min(score_cutoff);
            let cutoff_distance = maximum - score_cutoff;
            let hint_distance = maximum - score_hint;
            let dist = self._distance(s2_iter, len2, cutoff_distance, hint_distance);
            let sim = maximum - dist;
            if sim >= score_cutoff {
                sim
            } else {
                0 as $res_type
            }
        }
    };
}

macro_rules! build_cached_similarity_metric_funcs {
    ($impl_type:tt, $res_type:tt, $worst_similarity:expr, $worst_distance:expr) => {
        build_cached_normalized_metric_funcs!(
            $impl_type,
            $res_type,
            $worst_similarity,
            $worst_distance
        );

        #[allow(dead_code)]
        pub fn distance<Iter2, Elem2>(
            &self,
            s2: Iter2,
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>,
        ) -> $res_type
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let len2 = s2_iter.clone().count();
            self._distance(
                s2_iter,
                len2,
                score_cutoff.unwrap_or($worst_distance),
                score_hint.unwrap_or($worst_distance),
            )
        }

        #[allow(dead_code)]
        pub fn similarity<Iter2, Elem2>(
            &self,
            s2: Iter2,
            score_cutoff: Option<$res_type>,
            score_hint: Option<$res_type>,
        ) -> $res_type
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let len2 = s2_iter.clone().count();
            self._similarity(
                s2_iter,
                len2,
                score_cutoff.unwrap_or($worst_similarity),
                score_hint.unwrap_or($worst_similarity),
            )
        }

        pub(crate) fn _distance<Iter2, Elem2>(
            &self,
            s2: Iter2,
            len2: usize,
            score_cutoff: $res_type,
            score_hint: $res_type,
        ) -> $res_type
        where
            Iter2: IntoIterator<Item = Elem2>,
            Iter2::IntoIter: Clone,
            Elem1: PartialEq<Elem2> + HashableChar + Copy,
            Elem2: PartialEq<Elem1> + HashableChar + Copy,
            <Iter2 as IntoIterator>::IntoIter: DoubleEndedIterator,
        {
            let s2_iter = s2.into_iter();
            let maximum = self.maximum(s2_iter.clone(), len2);

            let cutoff_similarity = if maximum >= score_cutoff {
                maximum - score_cutoff
            } else {
                $worst_similarity as $res_type
            };
            let hint_similarity = if maximum >= score_hint {
                maximum - score_hint
            } else {
                $worst_similarity as $res_type
            };

            let sim = self._similarity(s2_iter, len2, cutoff_similarity, hint_similarity);
            let dist = maximum - sim;

            if dist <= score_cutoff {
                dist
            } else {
                less_than_score_cutoff_similarity!(score_cutoff, $res_type)
            }
        }
    };
}

pub(crate) use build_cached_distance_metric_funcs;
pub(crate) use build_cached_normalized_metric_funcs;
pub(crate) use build_cached_similarity_metric_funcs;
