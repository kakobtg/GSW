pub mod gadget;

use nalgebra::DMatrix;
use rand::Rng;
use crate::params::gsw::*;
use crate::FheError;
use gadget::{gadget_matrix, g_inverse};

#[derive(Clone, Debug)]
pub struct BitCiphertext(pub DMatrix<i64>);

#[derive(Clone, Debug)]
pub struct IntCiphertext(pub DMatrix<i64>);

pub fn keygen(rng: &mut impl Rng) -> (DMatrix<i64>, DMatrix<i64>) {
    let s = DMatrix::from_fn(N_LWE, 1, |_, _| rng.gen_range(0..=1_i64));
    let mut v = DMatrix::zeros(N_LWE + 1, 1);
    for i in 0..N_LWE {
        v[(i, 0)] = (-s[(i, 0)]).rem_euclid(Q);
    }
    v[(N_LWE, 0)] = 1;
    let b_mat = DMatrix::from_fn(N_LWE, M, |_, _| rng.gen_range(0..Q));
    let e = DMatrix::from_fn(1, M, |_, _| rng.gen_range(-1..=1));
    let b = (s.transpose() * &b_mat + e).map(|x| x.rem_euclid(Q));
    let mut a = DMatrix::zeros(N_LWE + 1, M);
    for i in 0..N_LWE {
        for j in 0..M { a[(i, j)] = b_mat[(i, j)]; }
    }
    for j in 0..M { a[(N_LWE, j)] = b[(0, j)]; }
    (v, a)
}

pub fn encrypt(a: &DMatrix<i64>, m: i64, rng: &mut impl Rng) -> Result<BitCiphertext, FheError> {
    if m != 0 && m != 1 {
        return Err(FheError::InvalidMessage { value: m, expected: "0 or 1" });
    }
    let r = DMatrix::from_fn(M, N_COLS, |_, _| rng.gen_range(0..=1));
    let g = gadget_matrix();
    let mut c = a * r + m * g;
    c.apply(|x| *x = x.rem_euclid(Q));
    Ok(BitCiphertext(c))
}

pub fn encrypt_integer(a: &DMatrix<i64>, m: i64, rng: &mut impl Rng) -> Result<IntCiphertext, FheError> {
    let max_val = 1 << INT_SCALE_SHIFT;
    if m.abs() >= max_val {
        return Err(FheError::ValueOutOfRange { value: m, max_abs: max_val });
    }
    let r = DMatrix::from_fn(M, N_COLS, |_, _| rng.gen_range(0..=1));
    let g = gadget_matrix();
    let mut c = a * r + m * g;
    c.apply(|x| *x = x.rem_euclid(Q));
    Ok(IntCiphertext(c))
}

pub fn decrypt(v: &DMatrix<i64>, c: &BitCiphertext) -> i64 {
    let vt_c = v.transpose() * &c.0;
    let val = vt_c[(0, N_COLS - 1)].rem_euclid(Q);
    if std::cmp::min(val, Q - val) < (val - (Q / 2)).abs() { 0 } else { 1 }
}

pub fn decrypt_integer(v: &DMatrix<i64>, c: &IntCiphertext) -> i64 {
    let vt_c = v.transpose() * &c.0;
    let val = vt_c[(0, N_LWE * L + INT_SCALE_SHIFT)].rem_euclid(Q);
    let centered_val = if val > Q / 2 { val - Q } else { val };
    (centered_val as f64 / (1_i64 << INT_SCALE_SHIFT) as f64).round() as i64
}

pub fn homomorphic_add(c1: &BitCiphertext, c2: &BitCiphertext) -> BitCiphertext {
    let mut c_sum = &c1.0 + &c2.0;
    c_sum.apply(|x| *x = x.rem_euclid(Q));
    BitCiphertext(c_sum)
}

pub fn homomorphic_add_int(c1: &IntCiphertext, c2: &IntCiphertext) -> IntCiphertext {
    let mut c_sum = &c1.0 + &c2.0;
    c_sum.apply(|x| *x = x.rem_euclid(Q));
    IntCiphertext(c_sum)
}

pub fn homomorphic_mul(c1: &BitCiphertext, c2: &BitCiphertext) -> BitCiphertext {
    let mut c_prod = &c1.0 * g_inverse(&c2.0);
    c_prod.apply(|x| *x = x.rem_euclid(Q));
    BitCiphertext(c_prod)
}

pub fn homomorphic_mul_int(c1: &IntCiphertext, c2: &IntCiphertext) -> IntCiphertext {
    let mut c_prod = &c1.0 * g_inverse(&c2.0);
    c_prod.apply(|x| *x = x.rem_euclid(Q));
    IntCiphertext(c_prod)
}

pub fn encrypt_bit_vector(a: &DMatrix<i64>, bits: &[i64], rng: &mut impl Rng) -> Result<Vec<BitCiphertext>, FheError> {
    bits.iter().map(|&b| encrypt(a, b, rng)).collect()
}

pub fn decrypt_bit_vector(v: &DMatrix<i64>, ciphertexts: &[BitCiphertext]) -> Vec<i64> {
    ciphertexts.iter().map(|c| decrypt(v, c)).collect()
}

pub fn homomorphic_equal(a: &DMatrix<i64>, vec1: &[BitCiphertext], vec2: &[BitCiphertext], rng: &mut impl Rng) -> Result<BitCiphertext, FheError> {
    if vec1.len() != vec2.len() {
        return Err(FheError::VectorLengthMismatch { len1: vec1.len(), len2: vec2.len() });
    }
    let enc_one = encrypt(a, 1, rng)?;
    let mut xnor_bits: Vec<_> = vec1.iter().zip(vec2.iter()).map(|(b1, b2)| {
        homomorphic_add(&homomorphic_add(b1, b2), &enc_one)
    }).collect();
    while xnor_bits.len() > 1 {
        let mut next_level = Vec::new();
        for chunk in xnor_bits.chunks(2) {
            if chunk.len() == 2 {
                next_level.push(homomorphic_mul(&chunk[0], &chunk[1]));
            } else {
                next_level.push(chunk[0].clone());
            }
        }
        xnor_bits = next_level;
    }
    Ok(xnor_bits.pop().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;
    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn encrypt_decrypt_roundtrip(m in 0i64..=1) {
            let mut rng = thread_rng();
            let (v, a) = keygen(&mut rng);
            let c = encrypt(&a, m, &mut rng).unwrap();
            prop_assert_eq!(decrypt(&v, &c), m);
        }
        #[test]
        fn homomorphic_add_is_xor(m1 in 0i64..=1, m2 in 0i64..=1) {
            let mut rng = thread_rng();
            let (v, a) = keygen(&mut rng);
            let c = homomorphic_add(&encrypt(&a, m1, &mut rng).unwrap(), &encrypt(&a, m2, &mut rng).unwrap());
            prop_assert_eq!(decrypt(&v, &c), m1 ^ m2);
        }
    }
}