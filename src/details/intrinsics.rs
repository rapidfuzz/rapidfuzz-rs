pub(crate) fn ceil_div_usize(a: usize, divisor: usize) -> usize {
    a / divisor + (a % divisor != 0) as usize
}

/// shift right without undefined behavior for shifts > bit width
pub(crate) fn shr64(a: u64, shift: usize) -> u64 {
    if shift < 64 {
        a >> shift
    } else {
        0
    }
}

/// shift left without undefined behavior for shifts > bit width
pub(crate) fn shl64(a: u64, shift: usize) -> u64 {
    if shift < 64 {
        a << shift
    } else {
        0
    }
}

// this is still a nightly only api. Can be removed if it becomes stable
pub(crate) fn carrying_add(lhs: u64, rhs: u64, carry: bool) -> (u64, bool) {
    let (a, b) = lhs.overflowing_add(rhs);
    let (c, d) = a.overflowing_add(carry as u64);
    (c, b | d)
}
