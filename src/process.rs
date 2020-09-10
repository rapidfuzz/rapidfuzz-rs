use crate::fuzz;

pub fn extract(query: &str, choices: &[&str], score_cutoff: f64) -> Vec<(String, f64)> {
	let mut result: Vec<(String, f64)> = choices.iter()
		.map(|choice| ((*choice).to_string(), fuzz::ratio(&query, &choice, score_cutoff)))
		.collect();
	result.sort_by(|(_, val_a), (_, val_b)| val_a.partial_cmp(val_b).unwrap());
	result
}

#[allow(non_snake_case)] // using non snake case for consistency with
pub fn extractOne(query: &str, choices: &[&str], mut score_cutoff: f64) -> Option<(String, f64)> {
	let mut result_choice = None;

	for choice in choices {
		let res = fuzz::ratio(&query, choice, score_cutoff);

		if res > score_cutoff {
			score_cutoff = res;
			result_choice = Some(((*choice).to_string(), res));
		}
	}

	result_choice
}