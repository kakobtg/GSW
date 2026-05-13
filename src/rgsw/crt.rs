use crate::params::rgsw::*;
use super::poly::Poly;

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

pub fn encode_slots(slots: &[i64; N]) -> Poly {
    let mut roots = [0; N];
    for i in 0..N { roots[i] = pow_mod(OMEGA, 2 * (i as i64) + 1, T); }
    let n_inv = pow_mod(N as i64, T - 2, T);
    let mut poly = Poly::zero();
    for j in 0..N {
        let mut sum = 0;
        for i in 0..N { sum = (sum + slots[i] * pow_mod(pow_mod(roots[i], T - 2, T), j as i64, T)) % T; }
        let mut coef = (sum * n_inv) % T;
        if coef > T / 2 { coef -= T; }
        poly.coefs[j] = coef;
    }
    poly
}

pub fn decode_slots(poly: &Poly) -> [i64; N] {
    let mut slots = [0; N];
    for i in 0..N {
        let mut val = 0;
        for j in 0..N { val = (val + poly.coefs[j].rem_euclid(T) * pow_mod(pow_mod(OMEGA, 2 * (i as i64) + 1, T), j as i64, T)) % T; }
        slots[i] = if val > T / 2 { val - T } else { val };
    }
    slots
}