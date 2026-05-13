use GSW::gsw::{
    keygen, encrypt, encrypt_integer, decrypt, decrypt_integer, homomorphic_add, 
    homomorphic_add_int, homomorphic_mul, homomorphic_mul_int, encrypt_bit_vector, 
    homomorphic_equal
};
use GSW::rgsw::{
    keygen as rgsw_keygen, encrypt as rgsw_encrypt, decrypt as rgsw_decrypt,
    homomorphic_add as rgsw_add, encrypt_poly, decrypt_poly, encrypt_integer as rgsw_enc_int,
    decrypt_integer as rgsw_dec_int, Poly
};
use rand::thread_rng;

fn run_bootstrapping_demo() {
    use GSW::params::rgsw::{N, Q};
    use GSW::bootstrapping::{generate_bootstrapping_key, bootstrap};
    use rand::Rng;

    println!("\n=================================================");
    println!("=== 3. Bootstrapping (Blind Rotation) Demo    ===");
    println!("=================================================");
    let mut rng = thread_rng();
    
    println!("[+] Generating RGSW Keys (Server Bootstrapping Keys)...");
    let (v_rgsw, a_rgsw) = rgsw_keygen(&mut rng);

    println!("[+] Generating LWE Client Key (Input Key)...");
    let mut lwe_sk = vec![0; N];
    for i in 0..4 { lwe_sk[i] = rng.gen_range(0..=1); }

    println!("[+] Generating Bootstrapping Key (Encrypted LWE secret)...");
    let bk = generate_bootstrapping_key(&a_rgsw, &lwe_sk, &mut rng);
    
    let q_4 = Q / 4;
    let mut test_poly = Poly::zero();
    for i in 0..N {
        if i < N / 4 || i > 3 * N / 4 { test_poly.coefs[i] = 0; } 
        else { test_poly.coefs[i] = (Q - q_4).rem_euclid(Q); }
    }
    println!("[+] Encrypting Test Polynomial (Lookup Table)...");
    let acc_init = encrypt_poly(&a_rgsw, &test_poly, &mut rng);

    let m = 1;
    println!("[+] Encrypting Message '{}' into LWE Ciphertext with massive noise...", m);
    let mut a = vec![0; N];
    let mut dot = 0i64;
    for i in 0..N {
        a[i] = rng.gen_range(0..Q);
        dot = ((dot as i128 + (a[i] as i128 * lwe_sk[i] as i128)) % Q as i128) as i64;
    }
    let b = (dot + rng.gen_range(-50..=50) + m * q_4).rem_euclid(Q);

    println!("[+] Scaling LWE Ciphertext from Z_Q down to Z_2N (Modulus Switching)...");
    let mut a_scaled = vec![0; N];
    for i in 0..N { a_scaled[i] = ((a[i] as f64 * (2.0 * N as f64) / (Q as f64)).round() as i64).rem_euclid(2 * N as i64); }
    let b_scaled = ((b as f64 * (2.0 * N as f64) / (Q as f64)).round() as i64).rem_euclid(2 * N as i64);

    println!("[+] Running Blind Rotation (Bootstrapping Circuit)...");
    let (ext_a, ext_b) = bootstrap(&acc_init, &bk, &a_scaled, b_scaled);
    
    println!("[+] Bootstrapping Complete! Decrypting Key-Switched Ciphertext...");
    let v0 = &v_rgsw.get(0, 0).coefs;
    let mut ext_sk = vec![0; N];
    ext_sk[0] = (-v0[0]).rem_euclid(Q);
    for i in 1..N { ext_sk[i] = v0[N - i].rem_euclid(Q); }
    let mut ext_dot = 0i64;
    for i in 0..N { ext_dot = ((ext_dot as i128 + (ext_a[i] as i128 * ext_sk[i] as i128)) % Q as i128) as i64; }
    
    let phase = (ext_b - ext_dot).rem_euclid(Q);
    let decrypted_m = if std::cmp::min(phase, Q - phase) < (phase - q_4).abs() { 0 } else { 1 };
    
    println!("    -> Extracted Message: {}", decrypted_m);
    if decrypted_m == m {
        println!("    -> Success! Noise cleared and message was perfectly recovered.");
    } else {
        println!("    -> Bootstrapping failed.");
    }
}

fn run_attack_demo() {
    use GSW::attack::export_lattice_basis_for_sage;
    println!("\n=================================================");
    println!("=== 4. Lattice Attack (Kannan's Embedding)    ===");
    println!("=================================================");
    let mut rng = thread_rng();
    let (_, a) = keygen(&mut rng);
    println!("[+] Generating fresh LWE Public Key...");
    println!("[+] Exporting LWE Public Key matrix to SageMath script format...");
    export_lattice_basis_for_sage(&a);
}

fn run_complex_circuit_with_bootstrapping_demo() {
    use GSW::params::rgsw::{N, Q, L, M};
    use GSW::rgsw::{Poly, PolyMatrix, encrypt as rgsw_encrypt, homomorphic_mul as rgsw_mul};
    use GSW::bootstrapping::{generate_bootstrapping_key, bootstrap};
    use rand::Rng;

    println!("\n=================================================");
    println!("=== 5. RGSW Complex Circuit & Bootstrapping   ===");
    println!("=================================================");
    let mut rng = thread_rng();

    println!("[+] Generating low-noise RGSW Keys for stable bootstrapping...");
    let mut s = Poly::zero();
    for i in 0..4 { s.coefs[i] = rng.gen_range(0..=1); } // 4 active bits for stability in toy Z_2N space
    let mut a_rgsw = PolyMatrix::zeros(M, 2);
    for i in 0..M {
        *a_rgsw.get_mut(i, 0) = Poly::rand_mod_q(&mut rng);
        *a_rgsw.get_mut(i, 1) = a_rgsw.get(i, 0).mul(&s).add(&Poly::rand_noise(&mut rng));
    }

    println!("[+] Encrypting inputs: m1=1, m2=1, m3=1...");
    let c1 = rgsw_encrypt(&a_rgsw, 1, &mut rng).unwrap();
    let c2 = rgsw_encrypt(&a_rgsw, 1, &mut rng).unwrap();
    let c3 = rgsw_encrypt(&a_rgsw, 1, &mut rng).unwrap();

    println!("[+] Evaluating Circuit: res = (m1 AND m2) AND m3");
    println!("    -> This depth-2 multiplication cascades error exponentially!");
    let c_and1 = rgsw_mul(&c1, &c2);
    let c_res = rgsw_mul(&c_and1, &c3);

    println!("[+] Extracting LWE Ciphertext from RGSW matrix...");
    // RGSW Column 2L-2 conveniently encrypts `m * Q/4` due to the Gadget Matrix structure.
    let target_col = 2 * L - 2; 
    let poly_a = c_res.0.get(0, target_col);
    let mut a_noisy = vec![0; N];
    for i in 0..N { a_noisy[i] = poly_a.coefs[i]; }
    let b_noisy = c_res.0.get(1, target_col).coefs[0];

    // Derive the exact LWE secret key from the underlying RGSW polynomial ring 
    let mut ext_sk = vec![0; N];
    ext_sk[0] = s.coefs[0];
    for i in 1..N { ext_sk[i] = (-s.coefs[N - i]).rem_euclid(Q); }

    let mut dot_noisy = 0i64;
    for i in 0..N { dot_noisy = ((dot_noisy as i128 + (a_noisy[i] as i128 * ext_sk[i] as i128)) % Q as i128) as i64; }
    let phase_noisy = (b_noisy - dot_noisy).rem_euclid(Q);
    let q_4 = Q / 4;
    
    let dist_noisy = (phase_noisy - q_4).rem_euclid(Q);
    let error_noisy = std::cmp::min(dist_noisy, Q - dist_noisy);
    println!("    -> Highly Noisy Phase: {} (Target: {})", phase_noisy, q_4);
    println!("    -> Accumulated Error : {}", error_noisy);

    println!("[+] Bootstrapping the damaged ciphertext to clear the noise...");
    
    // To bootstrap, we need the LWE ciphertext encrypted under a strictly binary key (s.coefs).
    // We can rewrite the extracted LWE `a` vector to naturally map to `s.coefs` instead of `ext_sk`.
    let mut a_prime = vec![0; N];
    a_prime[0] = a_noisy[0];
    for i in 1..N { a_prime[i] = (-a_noisy[N - i]).rem_euclid(Q); }

    let bk = generate_bootstrapping_key(&a_rgsw, &s.coefs, &mut rng);
    
    let mut test_poly = Poly::zero();
    for i in 0..N {
        if i < N / 4 || i > 3 * N / 4 { test_poly.coefs[i] = 0; } 
        else { test_poly.coefs[i] = (Q - q_4).rem_euclid(Q); }
    }
    let acc_init = GSW::rgsw::encrypt_poly(&a_rgsw, &test_poly, &mut rng);

    // Scale the a_prime vector for the blind rotation
    let mut a_scaled = vec![0; N];
    for i in 0..N { a_scaled[i] = ((a_prime[i] as f64 * (2.0 * N as f64) / (Q as f64)).round() as i64).rem_euclid(2 * N as i64); }
    let b_scaled = ((b_noisy as f64 * (2.0 * N as f64) / (Q as f64)).round() as i64).rem_euclid(2 * N as i64);

    let (a_clean, b_clean) = bootstrap(&acc_init, &bk, &a_scaled, b_scaled);

    // The bootstrap output is extracted natively from an RGSW matrix, meaning it is encrypted under ext_sk
    let mut dot_clean = 0i64;
    for i in 0..N { dot_clean = ((dot_clean as i128 + (a_clean[i] as i128 * ext_sk[i] as i128)) % Q as i128) as i64; }
    let phase_clean = (b_clean - dot_clean).rem_euclid(Q);
    
    let dist_clean = (phase_clean - q_4).rem_euclid(Q);
    let error_clean = std::cmp::min(dist_clean, Q - dist_clean);
    println!("    -> Bootstrapped Phase: {} (Target: {})", phase_clean, q_4);
    println!("    -> Residual Error    : {}", error_clean);

    if error_clean < error_noisy {
        println!("[+] Success! Noise was drastically reduced through Bootstrapping.");
    } else {
        println!("[-] Note: Due to toy parameters, noise might not perfectly reduce.");
    }
}

fn run_rgsw_demo() {
    let mut rng = thread_rng();
    let (v, a) = rgsw_keygen(&mut rng);
    let (m1, m3) = (1, 0);
    let (c1, c3) = (rgsw_encrypt(&a, m1, &mut rng).unwrap(), rgsw_encrypt(&a, m3, &mut rng).unwrap());

    println!("RGSW Fresh ciphertext decryption: m1={}, m3={}", rgsw_decrypt(&v, &c1), rgsw_decrypt(&v, &c3));
    println!("RGSW Homomorphic Add (m1+m3): Expected={}, Got={}", m1 ^ m3, rgsw_decrypt(&v, &rgsw_add(&c1, &c3)));

    println!("\n--- RGSW Polynomial Packing ---");
    let mut m_poly = Poly::zero();
    m_poly.coefs[0] = 1;
    m_poly.coefs[1] = 0;
    m_poly.coefs[2] = 1;
    
    let c_poly = encrypt_poly(&a, &m_poly, &mut rng);
    let dec_poly = decrypt_poly(&v, &c_poly);
    
    println!("Original packed bits: [1, 0, 1, 0, 0, 0, ...]");
    println!("Decrypted packed bits: [{}, {}, {}, {}, {}, {}, ...]", 
        dec_poly.coefs[0], dec_poly.coefs[1], dec_poly.coefs[2], 
        dec_poly.coefs[3], dec_poly.coefs[4], dec_poly.coefs[5]);

    println!("\n--- RGSW Scalar Small Integer Arithmetic ---");
    let int1 = 13;
    let int2 = 14;
    let c_int1 = rgsw_enc_int(&a, int1, &mut rng).unwrap();
    let c_int2 = rgsw_enc_int(&a, int2, &mut rng).unwrap();
    println!("RGSW Homomorphic Add ({} + {}): Got={}", int1, int2, rgsw_dec_int(&v, &rgsw_add(&c_int1, &c_int2)));
}

fn main() {
    println!("=================================================");
    println!("=== 1. GSW Homomorphic Encryption Test        ===");
    println!("=================================================");
    let mut rng = thread_rng();
    let (v, a) = keygen(&mut rng);

    let m1 = 1;
    let m2 = 1;
    let m3 = 0;

    println!("Encrypting messages: m1 = {}, m2 = {}, m3 = {}", m1, m2, m3);
    let c1 = encrypt(&a, m1, &mut rng).unwrap();
    let c2 = encrypt(&a, m2, &mut rng).unwrap();
    let c3 = encrypt(&a, m3, &mut rng).unwrap();

    assert_eq!(decrypt(&v, &c1), m1);
    assert_eq!(decrypt(&v, &c3), m3);
    println!("-> Fresh ciphertext decryption verified!");

    let c_add = homomorphic_add(&c1, &c3);
    println!("Homomorphic Addition (m1 + m3): Expected = {}, Got = {}", m1 ^ m3, decrypt(&v, &c_add));

    let c_mul_1 = homomorphic_mul(&c1, &c2);
    println!("Homomorphic Multiplication (m1 * m2): Expected = {}, Got = {}", m1 & m2, decrypt(&v, &c_mul_1));

    println!("\n=== Small Integer Arithmetic ===");
    let int1 = 13;
    let int2 = 14;
    println!("Encrypting integers: int1 = {}, int2 = {}", int1, int2);

    let c_int1 = encrypt_integer(&a, int1, &mut rng).unwrap();
    let c_int2 = encrypt_integer(&a, int2, &mut rng).unwrap();

    assert_eq!(decrypt_integer(&v, &c_int1), int1);
    println!("-> Fresh integer decryption verified!");

    let c_int_add = homomorphic_add_int(&c_int1, &c_int2);
    println!("Homomorphic Addition (int1 + int2): Expected = {}, Got = {}", int1 + int2, decrypt_integer(&v, &c_int_add));

    let c_int_mul = homomorphic_mul_int(&c_int1, &c_int2);
    println!("Homomorphic Multiplication (int1 * int2): Expected = {}, Got = {}", int1 * int2, decrypt_integer(&v, &c_int_mul));

    println!("\n=== Bit-Vector Evaluation ===");
    let vec1 = vec![1, 0, 1, 1];
    let vec2 = vec![1, 1, 0, 1];

    let c_vec1 = encrypt_bit_vector(&a, &vec1, &mut rng).unwrap();
    let c_vec2 = encrypt_bit_vector(&a, &vec2, &mut rng).unwrap();

    assert_eq!(decrypt(&v, &homomorphic_equal(&a, &c_vec1, &c_vec1, &mut rng).unwrap()), 1);
    assert_eq!(decrypt(&v, &homomorphic_equal(&a, &c_vec1, &c_vec2, &mut rng).unwrap()), 0);

    println!("\n=================================================");
    println!("=== 2. Ring-GSW (RGSW) Demo                   ===");
    println!("=================================================");
    run_rgsw_demo();

    run_bootstrapping_demo();
    
    run_attack_demo();

    run_complex_circuit_with_bootstrapping_demo();
}