pub mod gsw {
    pub const N_LWE: usize = 4;
    pub const Q: i64 = 1 << 24;
    pub const L: usize = 24;
    pub const N_COLS: usize = (N_LWE + 1) * L;
    pub const M: usize = 25;
    pub const INT_SCALE_SHIFT: usize = 14;
}

pub mod rgsw {
    pub const N: usize = 16;
    pub const Q: i64 = 1 << 35;
    pub const L: usize = 35;
    pub const M: usize = 4;
    pub const INT_SCALE_SHIFT: usize = 17;
    pub const T: i64 = 97;
    pub const OMEGA: i64 = 28;
}

#[cfg(test)]
mod param_checks {
    use super::rgsw::*;

    fn pow_mod(mut base: i64, mut exp: i64, modulus: i64) -> i64 {
        let mut res = 1;
        base = base.rem_euclid(modulus);
        while exp > 0 {
            if exp % 2 == 1 { res = (res * base) % modulus; }
            base = (base * base) % modulus;
            exp /= 2;
        }
        res
    }

    #[test]
    fn test_crt_params() {
        assert_eq!(pow_mod(OMEGA, 2 * N as i64, T), 1, "omega^2N must be 1 mod T");
        assert_ne!(pow_mod(OMEGA, N as i64, T), 1, "omega^N must not be 1 mod T");
    }
}