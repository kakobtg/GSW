use rand::Rng;
use std::cmp;

// --- Ring-GSW Cryptosystem Parameters ---
const N: usize = 16;                 // Polynomial degree (must be power of 2 for X^N+1)
const Q: i64 = 1 << 20;              // Modulus (q = 2^20)
const L: usize = 20;                 // log2(Q)
const M: usize = 4;                  // Number of RLWE samples in the public key

/// A polynomial in the ring Z_q[X] / (X^N + 1)
#[derive(Clone, Debug)]
pub struct Poly {
    pub coefs: [i64; N],
}

impl Poly {
    pub fn zero() -> Self { Poly { coefs: [0; N] } }
    
    pub fn from_scalar(s: i64) -> Self {
        let mut p = Poly::zero();
        p.coefs[0] = s.rem_euclid(Q);
        p
    }
    
    pub fn rand_mod_q() -> Self {
        let mut rng = rand::thread_rng();
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = rng.gen_range(0..Q); }
        p
    }
    
    pub fn rand_binary() -> Self {
        let mut rng = rand::thread_rng();
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = rng.gen_range(0..=1); }
        p
    }
    
    pub fn rand_noise() -> Self {
        let mut rng = rand::thread_rng();
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = rng.gen_range(-1..=1); }
        p
    }
    
    pub fn add(&self, other: &Poly) -> Self {
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = (self.coefs[i] + other.coefs[i]).rem_euclid(Q); }
        p
    }
    
    /// Naive polynomial multiplication modulo (X^N + 1)
    pub fn mul(&self, other: &Poly) -> Self {
        let mut p = Poly::zero();
        for i in 0..N {
            for j in 0..N {
                let coef = (self.coefs[i] * other.coefs[j]) % Q;
                if i + j >= N {
                    // Wrap around: X^N = -1
                    p.coefs[i + j - N] = (p.coefs[i + j - N] - coef).rem_euclid(Q);
                } else {
                    p.coefs[i + j] = (p.coefs[i + j] + coef).rem_euclid(Q);
                }
            }
        }
        p
    }
    
    pub fn scalar_mul(&self, scalar: i64) -> Self {
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = (self.coefs[i] * scalar).rem_euclid(Q); }
        p
    }
}

/// A Matrix where every element is a Polynomial
#[derive(Clone)]
pub struct PolyMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Poly>,
}

impl PolyMatrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        PolyMatrix { rows, cols, data: vec![Poly::zero(); rows * cols] }
    }
    pub fn get(&self, r: usize, c: usize) -> &Poly { &self.data[r * self.cols + c] }
    pub fn get_mut(&mut self, r: usize, c: usize) -> &mut Poly { &mut self.data[r * self.cols + c] }
    
    pub fn add(&self, other: &PolyMatrix) -> Self {
        let mut res = PolyMatrix::zeros(self.rows, self.cols);
        for i in 0..self.data.len() { res.data[i] = self.data[i].add(&other.data[i]); }
        res
    }
    
    pub fn mul(&self, other: &PolyMatrix) -> Self {
        assert_eq!(self.cols, other.rows, "Matrix dimension mismatch");
        let mut res = PolyMatrix::zeros(self.rows, other.cols);
        for r in 0..self.rows {
            for c in 0..other.cols {
                let mut sum = Poly::zero();
                for k in 0..self.cols {
                    sum = sum.add(&self.get(r, k).mul(other.get(k, c)));
                }
                *res.get_mut(r, c) = sum;
            }
        }
        res
    }
    
    pub fn transpose(&self) -> Self {
        let mut res = PolyMatrix::zeros(self.cols, self.rows);
        for r in 0..self.rows {
            for c in 0..self.cols { *res.get_mut(c, r) = self.get(r, c).clone(); }
        }
        res
    }
}

// --- RGSW Operations ---

pub fn keygen() -> (PolyMatrix, PolyMatrix) {
    let s = Poly::rand_binary(); 
    let mut v = PolyMatrix::zeros(2, 1);
    *v.get_mut(0, 0) = s.scalar_mul(-1); // -s
    *v.get_mut(1, 0) = Poly::from_scalar(1); // 1

    let mut a = PolyMatrix::zeros(M, 2);
    for i in 0..M {
        let a_poly = Poly::rand_mod_q();
        let e_poly = Poly::rand_noise();
        let b_poly = a_poly.mul(&s).add(&e_poly); // b = a * s + e
        *a.get_mut(i, 0) = a_poly;
        *a.get_mut(i, 1) = b_poly;
    }
    (v, a)
}

pub fn encrypt(a: &PolyMatrix, m: i64) -> PolyMatrix {
    let mut r = PolyMatrix::zeros(M, 2 * L);
    for i in 0..M {
        for j in 0..2 * L { *r.get_mut(i, j) = Poly::rand_binary(); }
    }
    
    let mut g = PolyMatrix::zeros(2, 2 * L);
    for i in 0..2 {
        for j in 0..L { *g.get_mut(i, i * L + j) = Poly::from_scalar(1 << j); }
    }
    
    let mut m_g = PolyMatrix::zeros(2, 2 * L);
    for i in 0..2 {
        for j in 0..2 * L { *m_g.get_mut(i, j) = g.get(i, j).scalar_mul(m); }
    }

    a.transpose().mul(&r).add(&m_g)
}

pub fn decrypt(v: &PolyMatrix, c: &PolyMatrix) -> i64 {
    let val_poly = v.transpose().mul(c).get(0, 2 * L - 1).clone();
    let val = val_poly.coefs[0].rem_euclid(Q); // Extract the scalar component
    if cmp::min(val, Q - val) < (val - (Q / 2)).abs() { 0 } else { 1 }
}

pub fn g_inverse(c: &PolyMatrix) -> PolyMatrix {
    let mut z = PolyMatrix::zeros(2 * L, 2 * L);
    for col in 0..2 * L {
        for row in 0..2 {
            let poly = c.get(row, col);
            for l in 0..L {
                let mut z_poly = Poly::zero();
                for d in 0..N { z_poly.coefs[d] = (poly.coefs[d] >> l) & 1; } // Decompose each coeff
                *z.get_mut(row * L + l, col) = z_poly;
            }
        }
    }
    z
}

pub fn homomorphic_add(c1: &PolyMatrix, c2: &PolyMatrix) -> PolyMatrix { c1.add(c2) }
pub fn homomorphic_mul(c1: &PolyMatrix, c2: &PolyMatrix) -> PolyMatrix { c1.mul(&g_inverse(c2)) }

pub fn run_rgsw_demo() {
    let (v, a) = keygen();
    let (m1, m2, m3) = (1, 1, 0);
    let (c1, c2, c3) = (encrypt(&a, m1), encrypt(&a, m2), encrypt(&a, m3));

    println!("RGSW Fresh ciphertext decryption: m1={}, m3={}", decrypt(&v, &c1), decrypt(&v, &c3));
    println!("RGSW Homomorphic Add (m1+m3): Expected={}, Got={}", m1 ^ m3, decrypt(&v, &homomorphic_add(&c1, &c3)));
    println!("RGSW Homomorphic Mul (m1*m2): Expected={}, Got={}", m1 & m2, decrypt(&v, &homomorphic_mul(&c1, &c2)));
    println!("RGSW Homomorphic Mul (m1*m3): Expected={}, Got={}", m1 & m3, decrypt(&v, &homomorphic_mul(&c1, &c3)));
}