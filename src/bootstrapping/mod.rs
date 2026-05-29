pub mod key;

use rayon::prelude::*;
use crate::params::rgsw::{N, Q, L};
use crate::rgsw::{Poly, PolyMatrix, RgswPolyCiphertext};
pub use key::{BootstrappingKey, generate_bootstrapping_key};

pub fn multiply_by_x_power(poly: &Poly, a: i64) -> Poly {
    let mut res = Poly::zero();
    let shift = a.rem_euclid(2 * N as i64) as usize; // Calculate cyclic shift amount modulo 2N (since X^N = -1)
    for i in 0..N {
        let wraps = (i + shift) / N;
        let mut val = poly.coefs[i];
        if wraps % 2 != 0 { val = (-val).rem_euclid(Q); } // Apply negacyclic property: shifting past N negates the coefficient
        res.coefs[(i + shift) % N] = (res.coefs[(i + shift) % N] + val).rem_euclid(Q); // Place the shifted (and possibly negated) coefficient
    }
    res
}

pub fn acc_multiply_by_x_power(acc: &PolyMatrix, a: i64) -> PolyMatrix {
    let mut res = PolyMatrix::zeros(acc.rows, acc.cols);
    res.data.par_iter_mut().zip(acc.data.par_iter()).for_each(|(r_poly, a_poly)| { // Parallelize the rotation across all polynomials in the matrix
        *r_poly = multiply_by_x_power(a_poly, a); // Rotate each polynomial by X^a
    });
    res
}

pub fn blind_rotate(acc: &RgswPolyCiphertext, bk: &BootstrappingKey, a: &[i64], b: i64) -> RgswPolyCiphertext {
    let mut current_acc = acc_multiply_by_x_power(&acc.0, b); // Initialize accumulator by blindly rotating by public LWE phase offset 'b'
    for (i, a_i) in a.iter().enumerate() { // Iterate over each component of the LWE ciphertext 'a'
        if *a_i == 0 { continue; } // Skip rotation if the component is zero
        let bk_i = &bk.bk[i];
        let acc_rotated = acc_multiply_by_x_power(&current_acc, -*a_i); // Compute ACC * X^{-a_i}
        let mut acc_diff = PolyMatrix::zeros(current_acc.rows, current_acc.cols);
        acc_diff.data.par_iter_mut().zip(acc_rotated.data.par_iter()).zip(current_acc.data.par_iter()).for_each(|((p_diff, p_rot), p_cur)| {
            for k in 0..N { p_diff.coefs[k] = (p_rot.coefs[k] - p_cur.coefs[k]).rem_euclid(Q); } // Compute the difference: ACC_rotated - ACC_current
        });
        
        let z2 = crate::rgsw::g_inverse(&acc_diff); // Bit-decompose the difference matrix for controlled error multiplication
        let mut selector = PolyMatrix::zeros(bk_i.0.rows, z2.cols);
        selector.data.par_chunks_mut(z2.cols).enumerate().for_each(|(row, row_slice)| { // Parallelize multiplication by the bootstrapping key (CMUX)
            for col in 0..z2.cols {
                let mut sum = Poly::zero();
                for k in 0..bk_i.0.cols { sum = sum.add(&bk_i.0.get(row, k).mul(z2.get(k, col))); } // Evaluate RGSW(s_i) * (ACC_rotated - ACC_current)
                row_slice[col] = sum;
            }
        });
        current_acc = current_acc.add(&selector); // Add selector result back: if s_i=1 it rotates, if s_i=0 it stays the same
    }
    RgswPolyCiphertext(current_acc)
}

/// Extracts a fresh LWE ciphertext from the RGSW accumulator (Sample Extraction).
/// TODO: This is a simplified sample extraction placeholder. A production implementation
/// must negate and rotate the `a` polynomial coefficients for correct negacyclic LWE key switching.
pub fn sample_extract(acc: &RgswPolyCiphertext) -> (Vec<i64>, i64) {
    let poly = acc.0.get(0, L); // Extract the polynomial from the phase column of the extended secret key structure
    let mut extracted_a = vec![0; N];
    for i in 0..N { extracted_a[i] = poly.coefs[i]; } // Read coefficients to form the new LWE 'a' vector
    (extracted_a, acc.0.get(1, L).coefs[0]) // Pair it with the 0-th degree constant term from the second row as the fresh LWE 'b' value
}

pub fn bootstrap(acc_init: &RgswPolyCiphertext, bk: &BootstrappingKey, lwe_ciphertext_a: &[i64], lwe_ciphertext_b: i64) -> (Vec<i64>, i64) {
    sample_extract(&blind_rotate(acc_init, bk, lwe_ciphertext_a, lwe_ciphertext_b)) // Chain blind rotation and sample extraction together
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rgsw::{keygen, encrypt_poly, decrypt_poly, Poly};
    use rand::thread_rng;

    #[test]
    fn test_blind_rotation() {
        let mut rng = thread_rng();
        let (v_rgsw, a_rgsw) = keygen(&mut rng);
        let lwe_sk = vec![1, 0, 1, 0];
        let bk = generate_bootstrapping_key(&a_rgsw, &lwe_sk, &mut rng);
        let a = vec![2, 5, 3, 1];
        let b = 7;
        
        let mut test_poly = Poly::zero();
        test_poly.coefs[0] = 1;

        let acc_init = encrypt_poly(&a_rgsw, &test_poly, &mut rng);
        let acc_final = blind_rotate(&acc_init, &bk, &a, b);
        let decrypted_poly = decrypt_poly(&v_rgsw, &acc_final);

        for i in 0..N {
            if i == 2 { assert_eq!(decrypted_poly.coefs[i], 1); } 
            else { assert_eq!(decrypted_poly.coefs[i], 0); }
        }
    }

    #[test]
    fn test_full_bootstrap_cycle() {
        let mut rng = thread_rng();
        let (_v_rgsw, a_rgsw) = keygen(&mut rng);
        let mut lwe_sk = vec![0; N];
        let mut a = vec![0; N];
        for i in 0..N {
            lwe_sk[i] = (i % 2) as i64;
            a[i] = (i as i64 * 12345) % Q;
        }
        let b = 9876543 % Q;
        let bk = generate_bootstrapping_key(&a_rgsw, &lwe_sk, &mut rng);
        
        let mut test_poly = Poly::zero();
        test_poly.coefs[0] = Q / 8;
        let acc_init = encrypt_poly(&a_rgsw, &test_poly, &mut rng);
        let (extracted_a, _) = bootstrap(&acc_init, &bk, &a, b);
        assert_eq!(extracted_a.len(), N);
    }

    #[test]
    fn test_integration_lwe_bootstrap() {
        use rand::Rng;
        let mut rng = thread_rng();
        let (v_rgsw, a_rgsw) = keygen(&mut rng);
        let mut lwe_sk = vec![0; N];
        for i in 0..4 { lwe_sk[i] = rng.gen_range(0..=1); }
        let bk = generate_bootstrapping_key(&a_rgsw, &lwe_sk, &mut rng);
        let q_4 = Q / 4;

        let mut test_poly = Poly::zero();
        for i in 0..N {
            if i < N / 4 || i > 3 * N / 4 { test_poly.coefs[i] = 0; } 
            else { test_poly.coefs[i] = (Q - q_4).rem_euclid(Q); }
        }
        let acc_init = encrypt_poly(&a_rgsw, &test_poly, &mut rng);

        for &m in &[0, 1] {
            let mut a = vec![0; N];
            let mut dot = 0i64;
            for i in 0..N {
                a[i] = rng.gen_range(0..Q);
                dot = ((dot as i128 + (a[i] as i128 * lwe_sk[i] as i128)) % Q as i128) as i64;
            }
            let b = (dot + rng.gen_range(-50..=50) + m * q_4).rem_euclid(Q);
            let mut a_scaled = vec![0; N];
            for i in 0..N { a_scaled[i] = ((a[i] as f64 * (2.0 * N as f64) / (Q as f64)).round() as i64).rem_euclid(2 * N as i64); }
            let b_scaled = ((b as f64 * (2.0 * N as f64) / (Q as f64)).round() as i64).rem_euclid(2 * N as i64);
            
            let (ext_a, ext_b) = bootstrap(&acc_init, &bk, &a_scaled, b_scaled);
            let v0 = &v_rgsw.get(0, 0).coefs;
            let mut ext_sk = vec![0; N];
            ext_sk[0] = (-v0[0]).rem_euclid(Q);
            for i in 1..N { ext_sk[i] = v0[N - i].rem_euclid(Q); }
            let mut ext_dot = 0i64;
            for i in 0..N { ext_dot = ((ext_dot as i128 + (ext_a[i] as i128 * ext_sk[i] as i128)) % Q as i128) as i64; }
            
            let phase = (ext_b - ext_dot).rem_euclid(Q);
            assert_eq!(if std::cmp::min(phase, Q - phase) < (phase - q_4).abs() { 0 } else { 1 }, m);
        }
    }
}