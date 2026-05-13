pub mod poly;
pub mod crt;

use rand::Rng;
use rayon::prelude::*;
use crate::params::rgsw::*;
use crate::FheError;
pub use poly::{Poly, PolyMatrix};
pub use crt::{encode_slots, decode_slots};

#[derive(Clone, Debug)]
pub struct RgswBitCiphertext(pub PolyMatrix);

#[derive(Clone, Debug)]
pub struct RgswPolyCiphertext(pub PolyMatrix);

pub fn keygen(rng: &mut impl Rng) -> (PolyMatrix, PolyMatrix) {
    let s = Poly::rand_binary(rng); 
    let mut v = PolyMatrix::zeros(2, 1);
    *v.get_mut(0, 0) = s.scalar_mul(-1);
    *v.get_mut(1, 0) = Poly::from_scalar(1);
    let mut a = PolyMatrix::zeros(M, 2);
    for i in 0..M {
        *a.get_mut(i, 0) = Poly::rand_mod_q(rng);
        *a.get_mut(i, 1) = a.get(i, 0).mul(&s).add(&Poly::rand_noise(rng));
    }
    (v, a)
}

pub fn encrypt(a: &PolyMatrix, m: i64, rng: &mut impl Rng) -> Result<RgswBitCiphertext, FheError> {
    if m != 0 && m != 1 {
        return Err(FheError::InvalidMessage { value: m, expected: "0 or 1" });
    }
    Ok(RgswBitCiphertext(core_encrypt(a, &Poly::from_scalar(m), rng)))
}

pub fn encrypt_integer(a: &PolyMatrix, m: i64, rng: &mut impl Rng) -> Result<RgswBitCiphertext, FheError> {
    let max_val = 1 << INT_SCALE_SHIFT;
    if m.abs() >= max_val {
        return Err(FheError::ValueOutOfRange { value: m, max_abs: max_val });
    }
    Ok(RgswBitCiphertext(core_encrypt(a, &Poly::from_scalar(m), rng)))
}

pub fn encrypt_poly(a: &PolyMatrix, m: &Poly, rng: &mut impl Rng) -> RgswPolyCiphertext {
    RgswPolyCiphertext(core_encrypt(a, m, rng))
}

fn core_encrypt(a: &PolyMatrix, m: &Poly, rng: &mut impl Rng) -> PolyMatrix {
    let mut r = PolyMatrix::zeros(M, 2 * L);
    for i in 0..M {
        for j in 0..2 * L { *r.get_mut(i, j) = Poly::rand_binary(rng); }
    }
    let mut m_g = PolyMatrix::zeros(2, 2 * L);
    for i in 0..2 {
        for j in 0..L { *m_g.get_mut(i, i * L + j) = Poly::from_scalar(1 << j).mul(m); }
    }
    
    // Implement `mul` internally for core operations, using parallel chunks over output rows
    let a_t = a.transpose();
    let mut res = PolyMatrix::zeros(a_t.rows, r.cols);
    res.data.par_chunks_mut(r.cols).enumerate().for_each(|(row, row_slice)| {
        for col in 0..r.cols {
            let mut sum = Poly::zero();
            for k in 0..a_t.cols { sum = sum.add(&a_t.get(row, k).mul(r.get(k, col))); }
            row_slice[col] = sum;
        }
    });
    res.add(&m_g)
}

pub fn decrypt(v: &PolyMatrix, c: &RgswBitCiphertext) -> i64 {
    let val = core_decrypt(v, &c.0, 2 * L - 1).coefs[0].rem_euclid(Q);
    if std::cmp::min(val, Q - val) < (val - (Q / 2)).abs() { 0 } else { 1 }
}

pub fn decrypt_integer(v: &PolyMatrix, c: &RgswBitCiphertext) -> i64 {
    let val = core_decrypt(v, &c.0, L + INT_SCALE_SHIFT).coefs[0].rem_euclid(Q);
    let centered = if val > Q / 2 { val - Q } else { val };
    (centered as f64 / (1_i64 << INT_SCALE_SHIFT) as f64).round() as i64
}

pub fn decrypt_poly(v: &PolyMatrix, c: &RgswPolyCiphertext) -> Poly {
    let mut res = Poly::zero();
    for (i, &val) in core_decrypt(v, &c.0, 2 * L - 1).coefs.iter().enumerate() {
        let val = val.rem_euclid(Q);
        res.coefs[i] = if std::cmp::min(val, Q - val) < (val - (Q / 2)).abs() { 0 } else { 1 };
    }
    res
}

pub fn encrypt_integer_poly(a: &PolyMatrix, m: &Poly, rng: &mut impl Rng) -> RgswBitCiphertext {
    RgswBitCiphertext(core_encrypt(a, m, rng))
}

pub fn decrypt_integer_poly(v: &PolyMatrix, c: &RgswBitCiphertext) -> Poly {
    let mut res = Poly::zero();
    for (i, &val) in core_decrypt(v, &c.0, L + INT_SCALE_SHIFT).coefs.iter().enumerate() {
        let val = val.rem_euclid(Q);
        let centered = if val > Q / 2 { val - Q } else { val };
        res.coefs[i] = (centered as f64 / (1_i64 << INT_SCALE_SHIFT) as f64).round() as i64;
    }
    res
}

fn core_decrypt(v: &PolyMatrix, c: &PolyMatrix, target_col: usize) -> Poly {
    let v_t = v.transpose();
    let mut res = Poly::zero();
    for k in 0..v_t.cols { res = res.add(&v_t.get(0, k).mul(c.get(k, target_col))); }
    res
}

pub fn g_inverse(c: &PolyMatrix) -> PolyMatrix {
    let mut z = PolyMatrix::zeros(2 * L, 2 * L);
    for col in 0..2 * L {
        for row in 0..2 {
            let poly = c.get(row, col);
            for l in 0..L {
                let mut z_poly = Poly::zero();
                for d in 0..N { z_poly.coefs[d] = (poly.coefs[d] >> l) & 1; }
                *z.get_mut(row * L + l, col) = z_poly;
            }
        }
    }
    z
}

pub fn homomorphic_add(c1: &RgswBitCiphertext, c2: &RgswBitCiphertext) -> RgswBitCiphertext { RgswBitCiphertext(c1.0.add(&c2.0)) }

pub fn homomorphic_mul(c1: &RgswBitCiphertext, c2: &RgswBitCiphertext) -> RgswBitCiphertext { 
    RgswBitCiphertext(c1.0.mul(&g_inverse(&c2.0))) 
}