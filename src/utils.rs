pub fn default_process(string: &str) -> String {

	let s: String = string.chars()
	  .map(|ch| 
		if ch <= '/' ||
		  (ch >= ':' && ch <= '@') ||
		  (ch >= '[' && ch <= '`') ||
		  (ch >= '{' && ch <= 0x7F as char) /* DEL */
		{
          ' '
		} else {
			ch
		}
	  ).collect();

	s.trim().to_ascii_lowercase()
}