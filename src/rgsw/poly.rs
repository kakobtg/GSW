use rand::Rng;
use rayon::prelude::*;
use crate::params::rgsw::*;

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
    
    pub fn rand_mod_q(rng: &mut impl Rng) -> Self {
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = rng.gen_range(0..Q); }
        p
    }
    
    pub fn rand_binary(rng: &mut impl Rng) -> Self {
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = rng.gen_range(0..=1); }
        p
    }
    
    pub fn rand_noise(rng: &mut impl Rng) -> Self {
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = rng.gen_range(-1..=1); }
        p
    }
    
    pub fn add(&self, other: &Poly) -> Self {
        let mut p = Poly::zero();
        for i in 0..N { p.coefs[i] = (self.coefs[i] + other.coefs[i]).rem_euclid(Q); }
        p
    }
    
    pub fn mul(&self, other: &Poly) -> Self {
        let mut p = Poly::zero();
        for i in 0..N {
            for j in 0..N {
                let prod = ((self.coefs[i] as i128 * other.coefs[j] as i128) % (Q as i128)) as i64;
                let target = i + j;
                if target >= N {
                    p.coefs[target - N] = (p.coefs[target - N] - prod).rem_euclid(Q);
                } else {
                    p.coefs[target] = (p.coefs[target] + prod).rem_euclid(Q);
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

#[derive(Clone, Debug)]
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
        res.data
            .par_chunks_mut(other.cols)
            .enumerate()
            .for_each(|(r, row_slice)| {
                for c in 0..other.cols {
                    let mut sum = Poly::zero();
                    for k in 0..self.cols {
                        sum = sum.add(&self.get(r, k).mul(other.get(k, c)));
                    }
                    row_slice[c] = sum;
                }
            });
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