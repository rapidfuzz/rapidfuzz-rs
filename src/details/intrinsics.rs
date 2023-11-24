pub(crate) const fn ceil_div_usize(a: usize, divisor: usize) -> usize {
    a / divisor + (a % divisor != 0) as usize
}

/// shift right without undefined behavior for shifts > bit width
pub(crate) const fn shr64(a: u64, shift: usize) -> u64 {
    if shift < 64 {
        a >> shift
    } else {
        0
    }
}

/// shift left without undefined behavior for shifts > bit width
#[allow(dead_code)]
pub(crate) const fn shl64(a: u64, shift: usize) -> u64 {
    if shift < 64 {
        a << shift
    } else {
        0
    }
}

// this is still a nightly only api. Can be removed if it becomes stable
pub(crate) const fn carrying_add(lhs: u64, rhs: u64, carry: bool) -> (u64, bool) {
    let (a, b) = lhs.overflowing_add(rhs);
    let (c, d) = a.overflowing_add(carry as u64);
    (c, b | d)
}

pub(crate) const fn bit_mask_lsb_u64(n: usize) -> u64 {
    let mut mask = !0_u64;
    if n < 64 {
        mask = mask.wrapping_add(1_u64 << n);
    }
    mask
}
pub(crate) const fn blsi_u64(v: u64) -> u64 {
    v & v.wrapping_neg()
}

pub(crate) const fn blsr_u64(v: u64) -> u64 {
    v & v.wrapping_sub(1)
}
