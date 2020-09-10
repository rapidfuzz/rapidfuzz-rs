pub fn remove_common_affix<'a, 'b>(first_string: &'a str, second_string: &'b str) -> (&'a str, &'b str) {
	// remove common prefix and suffix (linear vs square runtime for levensthein)
	let prefix_len = first_string.char_indices()
    	.zip(second_string.char_indices())
		.take_while(|&(a_char, b_char)| a_char == b_char)
		.count();

	let s1 = &first_string[prefix_len..];
	let s2 = &second_string[prefix_len..];

	let suffix_len = s1.chars().rev()
    	.zip(s2.chars().rev())
		.take_while(|&(a_char, b_char)| a_char == b_char)
		.count();

	(
		&s1[..s1.chars().count() - suffix_len],
		&s2[..s2.chars().count() - suffix_len]
	)
}


pub fn result_cutoff(result: f64, score_cutoff: f64) -> f64 {
	if result >= score_cutoff { result } else { 0.0 }
}

pub fn intersection_count_sorted_vec<'a>(mut a: Vec<&'a str>, mut b: Vec<&'a str>)
  -> (Vec<&'a str>, Vec<&'a str>, Vec<&'a str>)
{
    let mut sorted_sect: Vec<&str> = vec![];
    let mut sorted_1to2: Vec<&str> = vec![];
    a.dedup();
    b.dedup();

    for current_a in a {
        match b.binary_search(&current_a) {
            Ok(index) => {
                b.remove(index);
                sorted_sect.push(current_a)
            },
            _ => sorted_1to2.push(current_a),
        }
    }
    (sorted_sect, sorted_1to2, b)
}

pub fn sorted_split(sentence: &str) -> Vec<&str> {
	let mut tokens: Vec<_> = sentence.split_whitespace().collect();
	tokens.sort_unstable();
	tokens
}


pub fn norm_distance(dist: usize, lensum: usize, score_cutoff: f64) -> f64 {
  let ratio = 100.0 - 100.0 * dist as f64 / lensum as f64;
  result_cutoff(ratio, score_cutoff)
}

pub fn joined_size(sentence: Vec<&str>) -> usize {
	if sentence.is_empty() { return 0; }

	// there is a whitespace between each word
	let mut result = sentence.len() - 1;
	for word in sentence {
	  result += word.chars().count();
	}
  
	result
}
