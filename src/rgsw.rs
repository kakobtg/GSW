use rand::Rng;
use std::cmp;

// --- Ring-GSW Cryptosystem Parameters ---
const N: usize = 16;                 // Polynomial degree (must be power of 2 for X^N+1)
const Q: i64 = 1 << 35;              // Modulus (q = 2^35) to tolerate massive CRT coefficient growth
const L: usize = 35;                 // log2(Q)
const M: usize = 4;                  // Number of RLWE samples in the public key
const INT_SCALE_SHIFT: usize = 17;   // Column offset to use for small integer decryption
const T: i64 = 97;                   // Plaintext modulus for CRT packing (prime where T = 1 mod 2N)

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
                let coef = ((self.coefs[i] as i128 * other.coefs[j] as i128).rem_euclid(Q as i128)) as i64;
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
        for i in 0..N { p.coefs[i] = ((self.coefs[i] as i128 * scalar as i128).rem_euclid(Q as i128)) as i64; }
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

/// Encrypts a full polynomial (packing N bits) into a single RGSW ciphertext
pub fn encrypt_poly(a: &PolyMatrix, m: &Poly) -> PolyMatrix {
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
        for j in 0..2 * L { *m_g.get_mut(i, j) = g.get(i, j).mul(m); }
    }

    a.transpose().mul(&r).add(&m_g)
}

/// Decrypts an RGSW ciphertext back into a packed polynomial
pub fn decrypt_poly(v: &PolyMatrix, c: &PolyMatrix) -> Poly {
    let val_poly = v.transpose().mul(c).get(0, 2 * L - 1).clone();
    let mut res = Poly::zero();
    for i in 0..N {
        let val = val_poly.coefs[i].rem_euclid(Q);
        res.coefs[i] = if cmp::min(val, Q - val) < (val - (Q / 2)).abs() { 0 } else { 1 };
    }
    res
}

/// Encrypts a scalar integer into an RGSW ciphertext
pub fn encrypt_integer(a: &PolyMatrix, m: i64) -> PolyMatrix {
    encrypt(a, m) // Reuses the exact same encryption logic
}

/// Decrypts an RGSW ciphertext back into a scalar integer
pub fn decrypt_integer(v: &PolyMatrix, c: &PolyMatrix) -> i64 {
    let val_poly = v.transpose().mul(c).get(0, L + INT_SCALE_SHIFT).clone();
    let val = val_poly.coefs[0].rem_euclid(Q);
    
    let centered_val = if val > Q / 2 { val - Q } else { val };
    let scale = (1_i64 << INT_SCALE_SHIFT) as f64;
    
    (centered_val as f64 / scale).round() as i64
}

/// Encrypts a polynomial of small integers into an RGSW ciphertext
pub fn encrypt_integer_poly(a: &PolyMatrix, m: &Poly) -> PolyMatrix {
    encrypt_poly(a, m) // Reuses the exact same encryption logic
}

/// Decrypts an RGSW ciphertext back into a packed polynomial of small integers
pub fn decrypt_integer_poly(v: &PolyMatrix, c: &PolyMatrix) -> Poly {
    let val_poly = v.transpose().mul(c).get(0, L + INT_SCALE_SHIFT).clone();
    let mut res = Poly::zero();
    for i in 0..N {
        let val = val_poly.coefs[i].rem_euclid(Q);
        let centered_val = if val > Q / 2 { val - Q } else { val };
        let scale = (1_i64 << INT_SCALE_SHIFT) as f64;
        res.coefs[i] = (centered_val as f64 / scale).round() as i64;
    }
    res
}

/// Modular exponentiation helper
fn pow_mod(base: i64, exp: i64, modulus: i64) -> i64 {
    let mut res = 1;
    let mut b = base.rem_euclid(modulus);
    let mut e = exp;
    while e > 0 {
        if e % 2 == 1 { res = (res * b) % modulus; }
        b = (b * b) % modulus;
        e /= 2;
    }
    res
}

/// CRT Encode: packs N slots modulo T into a polynomial over Z_T[X]/(X^N+1)
pub fn encode_slots(slots: &[i64; N]) -> Poly {
    let omega = 28; // Primitive 32nd root of unity mod 97
    let mut roots = [0; N];
    for i in 0..N { roots[i] = pow_mod(omega, 2 * (i as i64) + 1, T); }

    let n_inv = pow_mod(N as i64, T - 2, T);
    let mut poly = Poly::zero();
    for j in 0..N {
        let mut sum = 0;
        for i in 0..N {
            let r_inv = pow_mod(roots[i], T - 2, T);
            let term = pow_mod(r_inv, j as i64, T);
            sum = (sum + slots[i] * term) % T;
        }
        let mut coef = (sum * n_inv) % T;
        if coef > T / 2 { coef -= T; }
        poly.coefs[j] = coef;
    }
    poly
}

/// CRT Decode: unpacks a polynomial evaluated at the roots modulo T back into N slots
pub fn decode_slots(poly: &Poly) -> [i64; N] {
    let omega = 28; // Primitive 32nd root of unity mod 97
    let mut roots = [0; N];
    for i in 0..N { roots[i] = pow_mod(omega, 2 * (i as i64) + 1, T); }

    let mut slots = [0; N];
    for i in 0..N {
        let mut val = 0;
        for j in 0..N {
            let coef = poly.coefs[j].rem_euclid(T);
            let x_pow = pow_mod(roots[i], j as i64, T);
            val = (val + coef * x_pow) % T;
        }
        if val > T / 2 { val -= T; }
        slots[i] = val;
    }
    slots
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

    println!("\n--- RGSW Polynomial Packing ---");
    // Create a polynomial message packing N bits
    let mut m_poly = Poly::zero();
    m_poly.coefs[0] = 1;
    m_poly.coefs[1] = 0;
    m_poly.coefs[2] = 1;
    m_poly.coefs[3] = 1;
    // The remaining 12 coefficients remain 0
    
    let c_poly = encrypt_poly(&a, &m_poly);
    let dec_poly = decrypt_poly(&v, &c_poly);
    
    println!("Original packed bits: [1, 0, 1, 1, 0, 0, ...]");
    println!("Decrypted packed bits: [{}, {}, {}, {}, {}, {}, ...]", 
        dec_poly.coefs[0], dec_poly.coefs[1], dec_poly.coefs[2], 
        dec_poly.coefs[3], dec_poly.coefs[4], dec_poly.coefs[5]);

    println!("\n--- RGSW Scalar Small Integer Arithmetic ---");
    let int1 = 13;
    let int2 = 14;
    let c_int1 = encrypt_integer(&a, int1);
    let c_int2 = encrypt_integer(&a, int2);
    println!("RGSW Homomorphic Add ({} + {}): Got={}", int1, int2, decrypt_integer(&v, &homomorphic_add(&c_int1, &c_int2)));
    println!("RGSW Homomorphic Mul ({} * {}): Got={}", int1, int2, decrypt_integer(&v, &homomorphic_mul(&c_int1, &c_int2)));

    println!("\n--- RGSW Element-Wise Polynomial Packing (CRT Slots) ---");
    let mut slots1 = [0; N];
    let mut slots2 = [0; N];
    
    slots1[0] = 5; slots1[1] = 12; slots1[2] = -4; slots1[3] = 7;
    slots2[0] = 4; slots2[1] = 3;  slots2[2] = 2;  slots2[3] = -2;
    
    let m_poly1 = encode_slots(&slots1);
    let m_poly2 = encode_slots(&slots2);
    
    let c_poly1 = encrypt_integer_poly(&a, &m_poly1);
    let c_poly2 = encrypt_integer_poly(&a, &m_poly2);
    
    let dec_add = decrypt_integer_poly(&v, &homomorphic_add(&c_poly1, &c_poly2));
    let slots_add = decode_slots(&dec_add);
    println!("Homomorphic Add (slots1 + slots2): [{}, {}, {}, {}, ...]", slots_add[0], slots_add[1], slots_add[2], slots_add[3]);
    
    let dec_mul = decrypt_integer_poly(&v, &homomorphic_mul(&c_poly1, &c_poly2));
    let slots_mul = decode_slots(&dec_mul);
    println!("Homomorphic Mul (slots1 * slots2) [Element-Wise CRT]: [{}, {}, {}, {}, ...]", slots_mul[0], slots_mul[1], slots_mul[2], slots_mul[3]);
}