use crate::levenshtein;
use crate::details;

use difflib::sequencematcher::SequenceMatcher;


pub fn ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    let result = levenshtein::normalized_weighted_distance(s1, s2, score_cutoff / 100.0);
    result * 100.0
}

pub fn partial_ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    if s1.is_empty() || s2.is_empty() || score_cutoff > 100.0 {
        return 0.0;
    }

    let (longer, shorter) = if s1.len() >= s2.len() {
        (s1, s2)
    } else {
        (s2, s1)
    };

    let mut matcher = SequenceMatcher::new(shorter, longer);
    let blocks = matcher.get_matching_blocks();

    let mut max_ratio = 0.0;
    for block in blocks {
        let long_start = if block.second_start > block.first_start {
            block.second_start - block.first_start
        } else {
            0
        };

        let long_end = long_start + shorter.chars().count();
        let long_substr = &longer[long_start..long_end];

        let ls_ratio = levenshtein::normalized_weighted_distance(shorter, long_substr, score_cutoff);
    
        if ls_ratio > 0.995 {
            return 100.0;
        }
        
        if ls_ratio > max_ratio {
            max_ratio = ls_ratio;
        }
            
    }

    details::result_cutoff(max_ratio * 100.0, score_cutoff)
}


pub fn token_sort_ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    if score_cutoff > 100.0 {
        return 0.0;
    }

    let tokens_a = details::sorted_split(s1);
    let tokens_b = details::sorted_split(s2);

    ratio(&tokens_a.join(" "), &tokens_b.join(" "), score_cutoff)
}


pub fn partial_token_sort_ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    if score_cutoff > 100.0 {
        return 0.0;
    }
    
    let tokens_a = details::sorted_split(s1);
    let tokens_b = details::sorted_split(s2);

    partial_ratio(&tokens_a.join(" "), &tokens_b.join(" "), score_cutoff)
}


pub fn token_set_ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    if score_cutoff > 100.0 {
        return 0.0;
    }
    
    let tokens_a = details::sorted_split(s1);
    let tokens_b = details::sorted_split(s2);

    let (sect, diff_ab, diff_ba) =
      details::intersection_count_sorted_vec(tokens_a, tokens_b);
    
    // one sentence is part of the other one
    if !sect.is_empty() && (diff_ab.is_empty() || diff_ba.is_empty()) {
      return 100.0;
    }

    let diff_ab_joined = diff_ab.join(" ");
    let diff_ba_joined = diff_ba.join(" ");

    let ab_len = diff_ab_joined.chars().count();
    let ba_len = diff_ba_joined.chars().count();
    let sect_len = details::joined_size(sect);

    // string length sect+ab <-> sect and sect+ba <-> sect
    let sect_ab_len = sect_len + !!sect_len as usize + ab_len;
    let sect_ba_len = sect_len + !!sect_len as usize + ba_len;

    let dist = levenshtein::weighted_distance(&diff_ab_joined, &diff_ba_joined);
    let result = details::norm_distance(dist, sect_ab_len + sect_ba_len, score_cutoff);

    if sect_len == 0 {
      return result;
    }

    // levenshtein distance sect+ab <-> sect and sect+ba <-> sect
    // since only sect is similar in them the distance can be calculated based on
    // the length difference
    let sect_ab_dist = !!sect_len as usize + ab_len;
    let sect_ab_ratio = details::norm_distance(sect_ab_dist, sect_len + sect_ab_len, score_cutoff);

    let sect_ba_dist = !!sect_len as usize + ba_len;
    let sect_ba_ratio = details::norm_distance(sect_ba_dist, sect_len + sect_ba_len, score_cutoff);

	result.max(sect_ab_ratio).max(sect_ba_ratio)
}


pub fn partial_token_set_ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    if score_cutoff > 100.0 {
        return 0.0;
    }
    
    let mut tokens_a: Vec<_> = s1.split_whitespace().collect();
    tokens_a.sort_unstable();
    let mut tokens_b: Vec<_> = s2.split_whitespace().collect();
    tokens_b.sort_unstable();

    let (sect, diff_ab, diff_ba) =
      details::intersection_count_sorted_vec(tokens_a, tokens_b);

    // exit early when there is a common word in both sequences
	if !sect.is_empty() {
		100.0
	} else {
		partial_ratio(&diff_ab.join(" "), &diff_ba.join(" "), score_cutoff)
	}
}

pub fn token_ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    if score_cutoff > 100.0 {
        return 0.0;
    }
    
    let mut tokens_a: Vec<_> = s1.split_whitespace().collect();
    tokens_a.sort_unstable();
    let mut tokens_b: Vec<_> = s2.split_whitespace().collect();
	tokens_b.sort_unstable();

    let (sect, diff_ab, diff_ba) =
      details::intersection_count_sorted_vec(tokens_a.clone(), tokens_b.clone());


    if !sect.is_empty() && (diff_ab.is_empty() || diff_ba.is_empty()) {
      return 100.0;
    }

	let diff_ab_joined = diff_ab.join(" ");
    let diff_ba_joined = diff_ba.join(" ");
    
    let ab_len = diff_ab_joined.chars().count();
    let ba_len = diff_ba_joined.chars().count();
	let sect_len = details::joined_size(sect);

    let mut result = ratio(&tokens_a.join(" "), &tokens_b.join(" "), score_cutoff);

    // string length sect+ab <-> sect and sect+ba <-> sect
    let sect_ab_len = sect_len + !!sect_len + ab_len;
    let sect_ba_len = sect_len + !!sect_len + ba_len;

    //TODO should use a different score_cutoff
    let dist = levenshtein::weighted_distance(&diff_ab_joined, &diff_ba_joined);
    result = result.max(details::norm_distance(dist, 2 * sect_ba_len, score_cutoff));

    // exit early since the other ratios are 0
    if sect_len == 0 {
        return result;
    }

    // levenshtein distance sect+ab <-> sect and sect+ba <-> sect
    // since only sect is similar in them the distance can be calculated based on
    // the length difference
    let sect_ab_dist = !!sect_len + ab_len;
    let sect_ab_ratio = details::norm_distance(sect_ab_dist, sect_len + sect_ab_len, score_cutoff);

    let sect_ba_dist = !!sect_len + ba_len;
    let sect_ba_ratio = details::norm_distance(sect_ba_dist, sect_len + sect_ba_len, score_cutoff);

    return result.max(sect_ab_ratio).max(sect_ba_ratio)
}

pub fn partial_token_ratio(s1: &str, s2: &str, score_cutoff: f64) -> f64 {
    if score_cutoff > 100.0 {
        return 0.0;
    }
    
    let mut tokens_a: Vec<_> = s1.split_whitespace().collect();
    tokens_a.sort_unstable();
    let mut tokens_b: Vec<_> = s2.split_whitespace().collect();
    tokens_b.sort_unstable();

    let (sect, diff_ab, diff_ba) =
      details::intersection_count_sorted_vec(tokens_a.clone(), tokens_b.clone());

    // exit early when there is a common word in both sequences
	if !sect.is_empty() {
		return 100.0;
    }
    
    let result = partial_ratio(&tokens_a.join(" "), &tokens_b.join(" "), score_cutoff);
    // do not calculate the same partial_ratio twice
    if tokens_a.len() == diff_ab.len() && tokens_b.len() == diff_ba.len() {
        return result;
    }

    result.max(partial_ratio(&diff_ab.join(" "), &diff_ba.join(" "), score_cutoff.max(result)))
}

#[allow(non_snake_case)] // using non snake case for consistency with
pub fn WRatio(s1: &str, s2: &str, mut score_cutoff: f64) -> f64 {
    if score_cutoff > 100.0 {
        return 0.0;
    }

    let unbase_scale = 0.95;

    let len_a = s1.chars().count();
    let len_b = s2.chars().count();
    let len_ratio = if len_a > len_b {
        len_a as f64 / len_b as f64
    } else {
        len_b as f64 / len_a as f64
    };

    let mut sratio = ratio(s1, s2, score_cutoff);
    score_cutoff = score_cutoff.max(sratio + 0.00001);

    if len_ratio < 1.5 {
        return sratio.max(token_ratio(s1, s2, score_cutoff / unbase_scale) * unbase_scale);
    }

    let partial_scale = if len_ratio < 8.0 { 0.9 } else { 0.6 };

    score_cutoff /= partial_scale;
    sratio = sratio.max(partial_ratio(s1, s2, score_cutoff) * partial_scale);

    // increase the score_cutoff by a small step so it might be able to exit early
    score_cutoff = score_cutoff.max(sratio + 0.00001) / unbase_scale;

    sratio.max(partial_token_ratio(s1, s2, score_cutoff) * unbase_scale * partial_scale)
}
