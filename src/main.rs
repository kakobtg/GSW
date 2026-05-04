use nalgebra::DMatrix;
use rand::Rng;

pub mod rgsw;

// --- GSW Cryptosystem Parameters ---
const N_LWE: usize = 4;                 // LWE dimension (n)
const Q: i64 = 1 << 24;                 // Modulus (q = 2^24) to tolerate deeper multiplication noise
const L: usize = 24;                    // log2(Q)
const N_COLS: usize = (N_LWE + 1) * L;  // N = (n+1) * log2(Q)
const M: usize = 25;                    // Number of LWE samples in the public key
const INT_SCALE_SHIFT: usize = 14;      // Column offset to use for small integer decryption (2^14) allows up to int 511

/// Generates the secret key `v` and public key matrix `A`
pub fn keygen() -> (DMatrix<i64>, DMatrix<i64>) {
    let mut rng = rand::thread_rng();

    // 1. Generate LWE secret s
    let s = DMatrix::from_fn(N_LWE, 1, |_, _| rng.gen_range(0..=1_i64)); // Use a short binary secret

    // 2. Construct secret eigenvector v = [-s, 1]^T
    let mut v = DMatrix::zeros(N_LWE + 1, 1);
    for i in 0..N_LWE {
        v[(i, 0)] = (-s[(i, 0)]).rem_euclid(Q);
    }
    v[(N_LWE, 0)] = 1;

    // 3. Construct public key A = [B; b] where b = s^T * B + e
    let b_mat = DMatrix::from_fn(N_LWE, M, |_, _| rng.gen_range(0..Q));
    let e = DMatrix::from_fn(1, M, |_, _| rng.gen_range(-1..=1)); // Small LWE noise

    let b = (s.transpose() * &b_mat + e).map(|x| x.rem_euclid(Q));

    let mut a = DMatrix::zeros(N_LWE + 1, M);
    for i in 0..N_LWE {
        for j in 0..M {
            a[(i, j)] = b_mat[(i, j)];
        }
    }
    for j in 0..M {
        a[(N_LWE, j)] = b[(0, j)];
    }

    (v, a)
}

/// Constructs the Gadget matrix G (Powers of 2 layout for bit decomposition)
pub fn gadget_matrix() -> DMatrix<i64> {
    let mut g = DMatrix::zeros(N_LWE + 1, N_COLS);
    for i in 0..=N_LWE {
        for j in 0..L {
            g[(i, i * L + j)] = 1 << j;
        }
    }
    g
}

/// Encrypts a bit (0 or 1) into an (n+1) x N matrix: C = A * R + m * G
pub fn encrypt(a: &DMatrix<i64>, m: i64) -> DMatrix<i64> {
    assert!(m == 0 || m == 1, "Message must be 0 or 1");
    let mut rng = rand::thread_rng();
    
    // R is a small random matrix in {0, 1}
    let r = DMatrix::from_fn(M, N_COLS, |_, _| rng.gen_range(0..=1));
    let g = gadget_matrix();

    let mut c = a * r + m * g;
    c.apply(|x| *x = x.rem_euclid(Q));
    c
}

/// Encrypts a small integer into a ciphertext matrix
pub fn encrypt_integer(a: &DMatrix<i64>, m: i64) -> DMatrix<i64> {
    let mut rng = rand::thread_rng();
    
    // R is a small random matrix in {0, 1}
    let r = DMatrix::from_fn(M, N_COLS, |_, _| rng.gen_range(0..=1));
    let g = gadget_matrix();

    let mut c = a * r + m * g;
    c.apply(|x| *x = x.rem_euclid(Q));
    c
}

/// Decrypts a ciphertext matrix back into a small integer
pub fn decrypt_integer(v: &DMatrix<i64>, c: &DMatrix<i64>) -> i64 {
    let vt_c = v.transpose() * c;
    
    // Extract from the column associated with our chosen scaling factor (2^INT_SCALE_SHIFT)
    let val = vt_c[(0, N_LWE * L + INT_SCALE_SHIFT)].rem_euclid(Q);
    
    // Center the value around 0 to gracefully handle LWE noise
    let centered_val = if val > Q / 2 { val - Q } else { val };
    let scale = (1_i64 << INT_SCALE_SHIFT) as f64;
    
    (centered_val as f64 / scale).round() as i64
}

/// Decrypts the ciphertext matrix C back into a bit
pub fn decrypt(v: &DMatrix<i64>, c: &DMatrix<i64>) -> i64 {
    // Left-multiply by the secret vector: v^T * C
    let vt_c = v.transpose() * c;

    // We extract the message from the very last element of the vector.
    // In the gadget matrix G, the last element corresponds to 1 * 2^{L-1} (which is Q/2).
    // Therefore, the last element of v^T * C ≈ m * (Q/2).
    let val = vt_c[(0, N_COLS - 1)].rem_euclid(Q);

    // Check if the value is closer to 0 or Q/2
    let dist0 = std::cmp::min(val, Q - val);
    let dist_half = (val - (Q / 2)).abs();

    if dist0 < dist_half {
        0
    } else {
        1
    }
}

/// Homomorphic Addition: C_sum = C1 + C2
pub fn homomorphic_add(c1: &DMatrix<i64>, c2: &DMatrix<i64>) -> DMatrix<i64> {
    let mut c_sum = c1 + c2;
    c_sum.apply(|x| *x = x.rem_euclid(Q));
    c_sum
}

/// Applies the G_inverse / Flatten operation on Matrix C (Bit Decomposition)
pub fn g_inverse(c: &DMatrix<i64>) -> DMatrix<i64> {
    let mut z = DMatrix::zeros(N_COLS, N_COLS);
    
    // Decompose every element into its base-2 representation
    for col in 0..N_COLS {
        for row in 0..=N_LWE {
            let mut val = c[(row, col)].rem_euclid(Q);
            for j in 0..L {
                z[(row * L + j, col)] = val % 2;
                val /= 2;
            }
        }
    }
    z
}

/// Homomorphic Multiplication: C_prod = C1 * G_inverse(C2)
pub fn homomorphic_mul(c1: &DMatrix<i64>, c2: &DMatrix<i64>) -> DMatrix<i64> {
    let z2 = g_inverse(c2);
    let mut c_prod = c1 * z2;
    c_prod.apply(|x| *x = x.rem_euclid(Q));
    c_prod
}

/// Encrypts an array of bits into a vector of ciphertexts
pub fn encrypt_bit_vector(a: &DMatrix<i64>, bits: &[i64]) -> Vec<DMatrix<i64>> {
    bits.iter().map(|&b| encrypt(a, b)).collect()
}

/// Decrypts a vector of ciphertexts back into an array of bits
pub fn decrypt_bit_vector(v: &DMatrix<i64>, ciphertexts: &[DMatrix<i64>]) -> Vec<i64> {
    ciphertexts.iter().map(|c| decrypt(v, c)).collect()
}

/// Evaluates if two encrypted bit-vectors are equal. Returns an encrypted 1 if equal, 0 otherwise.
pub fn homomorphic_equal(a: &DMatrix<i64>, vec1: &[DMatrix<i64>], vec2: &[DMatrix<i64>]) -> DMatrix<i64> {
    assert_eq!(vec1.len(), vec2.len(), "Bit-vectors must be of the same length");
    
    let enc_one = encrypt(a, 1);

    // 1. Bitwise XNOR: NOT(A XOR B)
    let mut xnor_bits: Vec<_> = vec1.iter().zip(vec2.iter()).map(|(b1, b2)| {
        let xor_bit = homomorphic_add(b1, b2);
        homomorphic_add(&xor_bit, &enc_one) // NOT operation
    }).collect();

    // 2. AND all XNOR bits together using a balanced binary tree to control multiplicative depth
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

    xnor_bits.pop().unwrap()
}

/// Exports the public key into a lattice basis matrix formatted for a SageMath attack script.
/// This uses Kannan's Embedding to turn the LWE Bounded Distance Decoding problem into a Shortest Vector Problem.
pub fn export_lattice_basis_for_sage(a: &DMatrix<i64>) {
    println!("\n# === SAGE MATH LATTICE ATTACK SCRIPT ===");
    println!("# Copy and paste everything below this line into SageMath (or an online Sage cell) to run the attack:");
    println!("m = {}", M);
    println!("n = {}", N_LWE);
    println!("q = {}", Q);
    println!("basis = [");
    
    // 1. q * I_m rows
    for i in 0..M {
        let mut row = vec![0; M + N_LWE + 1];
        row[i] = Q;
        println!("    {:?},", row);
    }
    // 2. B^T and I_n rows
    for i in 0..N_LWE {
        let mut row = vec![0; M + N_LWE + 1];
        for j in 0..M { row[j] = a[(i, j)]; }
        row[M + i] = 1;
        println!("    {:?},", row);
    }
    // 3. -b^T row (Kannan's embedding)
    {
        let mut row = vec![0; M + N_LWE + 1];
        for j in 0..M { row[j] = -a[(N_LWE, j)]; }
        row[M + N_LWE] = 1;
        println!("    {:?}", row); 
    }
    println!("]");
    
    println!("L = Matrix(ZZ, basis)");
    println!("print(\"Running LLL / BKZ lattice reduction...\")");
    println!("L_red = L.LLL() # Or use L.BKZ(block_size=10) for stronger reduction");
    println!("for row in L_red:");
    println!("    if row[-1] == 1 or row[-1] == -1:");
    println!("        secret = row[m:m+n]");
    println!("        if row[-1] == -1: secret = [-x for x in secret]");
    println!("        print(\"\\n[+] LWE Broken! Recovered Secret s:\", secret)");
    println!("print(\"\\n[+] FINISHED\")");

}

/// Test Script mapping out evaluations
fn main() {
    println!("=== GSW Homomorphic Encryption Test ===");
    
    let (v, a) = keygen();

    let m1 = 1;
    let m2 = 1;
    let m3 = 0;

    println!("Encrypting messages: m1 = {}, m2 = {}, m3 = {}", m1, m2, m3);
    let c1 = encrypt(&a, m1);
    let c2 = encrypt(&a, m2);
    let c3 = encrypt(&a, m3);

    //println!("\nSample Ciphertext c1 (Dimensions: {} x {}):\n{}", c1.nrows(), c1.ncols(), c1);

    // Test Decryption on fresh ciphertexts
    assert_eq!(decrypt(&v, &c1), m1);
    assert_eq!(decrypt(&v, &c3), m3);
    println!("-> Fresh ciphertext decryption verified!");

    // Test Homomorphic Addition (Acts as XOR mod 2 because messages map to {0, Q/2})
    let c_add = homomorphic_add(&c1, &c3); // 1 + 0 = 1
    let dec_add = decrypt(&v, &c_add);
    println!("Homomorphic Addition (m1 + m3): Expected = {}, Got = {}", m1 ^ m3, dec_add);
    assert_eq!(dec_add, m1 ^ m3);

    // Test Homomorphic Multiplication (Acts as AND)
    let c_mul_1 = homomorphic_mul(&c1, &c2); // 1 * 1 = 1
    let dec_mul_1 = decrypt(&v, &c_mul_1);
    println!("Homomorphic Multiplication (m1 * m2): Expected = {}, Got = {}", m1 & m2, dec_mul_1);
    assert_eq!(dec_mul_1, m1 & m2);

    let c_mul_2 = homomorphic_mul(&c1, &c3); // 1 * 0 = 0
    let dec_mul_2 = decrypt(&v, &c_mul_2);
    println!("Homomorphic Multiplication (m1 * m3): Expected = {}, Got = {}", m1 & m3, dec_mul_2);
    assert_eq!(dec_mul_2, m1 & m3);

    println!("-> All homomorphic operations successfully evaluated!");

    println!("\n=== Small Integer Arithmetic ===");
    let int1 = 13;
    let int2 = 14;
    println!("Encrypting integers: int1 = {}, int2 = {}", int1, int2);

    let c_int1 = encrypt_integer(&a, int1);
    let c_int2 = encrypt_integer(&a, int2);

    assert_eq!(decrypt_integer(&v, &c_int1), int1);
    assert_eq!(decrypt_integer(&v, &c_int2), int2);
    println!("-> Fresh integer decryption verified!");

    let c_int_add = homomorphic_add(&c_int1, &c_int2);
    println!("Homomorphic Addition (int1 + int2): Expected = {}, Got = {}", int1 + int2, decrypt_integer(&v, &c_int_add));

    let c_int_mul = homomorphic_mul(&c_int1, &c_int2);
    println!("Homomorphic Multiplication (int1 * int2): Expected = {}, Got = {}", int1 * int2, decrypt_integer(&v, &c_int_mul));

    println!("\n=== Bit-Vector Evaluation ===");
    let vec1 = vec![1, 0, 1, 1];
    let vec2 = vec![1, 1, 0, 1];

    println!("Encrypting bit-vectors:");
    println!("vec1 = {:?}", vec1);
    println!("vec2 = {:?}", vec2);

    let c_vec1 = encrypt_bit_vector(&a, &vec1);
    let c_vec2 = encrypt_bit_vector(&a, &vec2);

    let dec_vec1 = decrypt_bit_vector(&v, &c_vec1);
    assert_eq!(dec_vec1, vec1);
    println!("-> Fresh bit-vector decryption verified!");

    // Test bitwise XOR (Addition)
    let c_vec_xor: Vec<_> = c_vec1.iter().zip(c_vec2.iter()).map(|(c1, c2)| homomorphic_add(c1, c2)).collect();
    let expected_xor: Vec<_> = vec1.iter().zip(vec2.iter()).map(|(b1, b2)| b1 ^ b2).collect();
    let dec_vec_xor = decrypt_bit_vector(&v, &c_vec_xor);
    println!("Homomorphic Bitwise XOR: Expected = {:?}, Got = {:?}", expected_xor, dec_vec_xor);
    assert_eq!(dec_vec_xor, expected_xor);

    // Test bitwise AND (Multiplication)
    let c_vec_and: Vec<_> = c_vec1.iter().zip(c_vec2.iter()).map(|(c1, c2)| homomorphic_mul(c1, c2)).collect();
    let expected_and: Vec<_> = vec1.iter().zip(vec2.iter()).map(|(b1, b2)| b1 & b2).collect();
    let dec_vec_and = decrypt_bit_vector(&v, &c_vec_and);
    println!("Homomorphic Bitwise AND: Expected = {:?}, Got = {:?}", expected_and, dec_vec_and);
    assert_eq!(dec_vec_and, expected_and);

    // Test Equality Check
    let c_eq1 = homomorphic_equal(&a, &c_vec1, &c_vec1);
    println!("Homomorphic Equality (vec1 == vec1): Expected = 1, Got = {}", decrypt(&v, &c_eq1));
    assert_eq!(decrypt(&v, &c_eq1), 1);

    let c_eq2 = homomorphic_equal(&a, &c_vec1, &c_vec2);
    println!("Homomorphic Equality (vec1 == vec2): Expected = 0, Got = {}", decrypt(&v, &c_eq2));
    assert_eq!(decrypt(&v, &c_eq2), 0);

    // Export the LWE problem for a SageMath Attack
    //export_lattice_basis_for_sage(&a);
    //println!("\nOriginal Secret s in Rust was: {:?}", v.view((0, 0), (N_LWE, 1)).iter().map(|x| (-x).rem_euclid(Q)).collect::<Vec<_>>());

    println!("\n=== Running Ring-GSW (RGSW) Demo ===");
    rgsw::run_rgsw_demo();
}
