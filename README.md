# GSW Homomorphic Encryption in Rust

A functional implementation of the Gentry-Sahai-Waters (GSW) homomorphic encryption scheme (GSW13) in Rust. This project demonstrates Key Generation, Encryption, Decryption, and Homomorphic operations (Addition and Multiplication) based on the Learning with Errors (LWE) problem. It also includes an implementation of **Ring-GSW (RGSW)** based on Ring-LWE.

## Prerequisites

- [Rust](https://www.rust-lang.org/tools/install) (includes Cargo)

The project relies on the following crates (defined in `Cargo.toml`):
- `nalgebra` for matrix operations.
- `rand` for random number generation.
- `rayon` for data-parallelism during bootstrapping.
- `proptest` for property-based cryptographic testing.

## Mathematical Foundation

As described in the GSW13 paper, the core of the scheme relies on the **approximate eigenvector method**. In this framework:
- A **ciphertext** $C$ is a matrix.
- The **secret key** $\vec{v}$ acts as an approximate eigenvector.
- The **message** $\mu$ is the eigenvalue.

The decryption relation is defined as:
$C \cdot \vec{v} = \mu \cdot \vec{v} + \vec{e}$
where $\vec{e}$ is a "small" error vector introduced for LWE security.

### Homomorphic Addition
Matrix addition corresponds to homomorphic addition of the underlying messages:
$(C_1 + C_2) \cdot \vec{v} = C_1 \cdot \vec{v} + C_2 \cdot \vec{v} = (\mu_1 + \mu_2) \cdot \vec{v} + (\vec{e}_1 + \vec{e}_2)$
The error grows additively, which is standard and easily tolerated.

### Homomorphic Multiplication
Matrix multiplication corresponds to homomorphic multiplication:
$C_1 \cdot C_2 \cdot \vec{v} = C_1 \cdot (\mu_2 \cdot \vec{v} + \vec{e}_2) = \mu_2 \cdot (C_1 \cdot \vec{v}) + C_1 \cdot \vec{e}_2 = \mu_1 \cdot \mu_2 \cdot \vec{v} + \mu_2 \cdot \vec{e}_1 + C_1 \cdot \vec{e}_2$

To prevent the term $C_1 \cdot \vec{e}_2$ from causing exponential error growth, the scheme ensures that the matrices used in multiplication have very small coefficients (e.g., $0$ or $1$) by utilizing the Gadget Matrix and Bit Decomposition.

## Ring-GSW (RGSW)

This project also features an RGSW implementation (`src/rgsw/mod.rs`) which transitions from standard LWE to Ring-LWE. 

### Mathematical Foundation of RGSW
Instead of operating on matrices of integers, RGSW operates on matrices of **polynomials** in the ring $R_q = Z_q[X]/(X^N + 1)$. 
- **Polynomial Matrices:** Ciphertexts are matrices where every element is a polynomial in $R_q$.
- **Secret Key:** The secret key is a vector of polynomials, typically $\vec{v} = (-s, 1)^T$, where $s$ is a binary or small polynomial.
- **Decryption Relation:** The approximate eigenvector method holds identically over the polynomial ring:
  $C \cdot \vec{v} = \mu \cdot \vec{v} + \vec{e}$
  where multiplication is polynomial multiplication modulo $(X^N + 1)$ and modulo $q$.

### Advantages and Coefficient Bit-Decomposition
- **Smaller Dimensions:** Because a single polynomial of degree $N$ effectively packs $N$ LWE samples into one object, matrix dimensions shrink from $(n+1) \times (n+1)\log_2(q)$ to just $2 \times 2\log_2(q)$.
- **Bit-Decomposition ($G^{-1}$):** The $G^{-1}$ flattening function operates over every individual coefficient of the polynomial independently. For a polynomial $p(X) = \sum p_i X^i$, the bit decomposition outputs polynomials where every coefficient is strictly $0$ or $1$, ensuring noise growth remains strictly bounded during homomorphic multiplication.

## Project Structure

The codebase is organized into a clean, idiomatic Rust library and binary structure:
- **`src/lib.rs`**: Core library entry point and error definitions.
- **`src/params.rs`**: Centralized cryptographic parameters for both GSW and RGSW.
- **`src/gsw/`**: Standard GSW implementation (LWE-based) including encryption, decryption, homomorphic ops, and the gadget matrix.
- **`src/rgsw/`**: Ring-GSW implementation (Ring-LWE-based) including polynomial packing and CRT slots.
- **`src/bootstrapping/`**: TFHE/FHEW Blind Rotation and Bootstrapping Key generation.
- **`src/attack.rs`**: Exports a given LWE public key matrix into a SageMath script for Lattice Reduction attacks.
- **`src/bin/main.rs`**: The comprehensive demo runner showcasing all FHE capabilities.

## Running the Demo and Tests

### 1. Run the Comprehensive Demo
To execute the full demonstration covering GSW operations, Ring-GSW polynomial packing, Bootstrapping, and the Lattice Attack export, run:

```bash
cargo run --bin main
```

This will print a step-by-step walkthrough to your console. It will also automatically generate a `lattice_attack.py` script in your directory, which you can run via SageMath:
```bash
sage lattice_attack.py
```

### 2. Run the Cryptographic Test Suite
The project uses standard Rust tests alongside **`proptest`** to mathematically verify homomorphic properties over thousands of random permutations.

To run all tests (including the bootstrapping cycles):
```bash
cargo test
```

To run a deeper cryptographic check (e.g., 1000 iterations for property tests):
```bash
PROPTEST_CASES=1000 cargo test
```

To view the statistics and timers during tests, disable output capture:
```bash
cargo test -- --nocapture
```

## Implemented Operations

- **KeyGen:** Generates a secret vector `v` and a public matrix `A` based on LWE.
- **Encrypt:** Encodes a message bit as an eigenvalue of a randomized matrix: `C = A * R + m * G`.
- **Decrypt:** Recovers the message by computing `v^T * C` and analyzing the extracted scalar.
- **Homomorphic Addition:** Simply adds two ciphertext matrices `C1 + C2` (acts as XOR).
- **Homomorphic Multiplication:** Multiplies ciphertexts using the bit-decomposition method `C1 * G_inverse(C2)` (acts as AND).

## Noise Growth and `G_inverse`

In standard LWE-based schemes, multiplying two ciphertexts multiplies their respective noise terms ($E_1 \cdot E_2$), causing an exponential blowup in noise and quickly destroying the underlying message.

The GSW scheme manages this using the **Gadget Matrix ($G$)** and the **Flattening/Bit-Decomposition function ($G^{-1}$)**:

1. **Gadget Matrix ($G$):** A matrix populated with powers of 2.
2. **Bit Decomposition ($G^{-1}$):** When applying $G^{-1}(C_2)$, the ciphertext matrix is decomposed into a binary matrix containing only `0`s and `1`s.

### Multiplication Noise Trace
When computing $C_{prod} = C_1 \cdot G^{-1}(C_2)$, the decryption relation maps to:

$v^T \cdot C_{prod} \approx [E_1 \cdot G^{-1}(C_2) + m_1 \cdot E_2] + (m_1 \cdot m_2) \cdot v^T \cdot G$

Because $G^{-1}(C_2)$ consists purely of binary values ($\{0, 1\}$), multiplying it by the noise $E_1$ acts as a simple addition of the noise terms. As a result, the noise grows **asymmetrically and linearly** (scaled at worst by the matrix dimension $N$), rather than exponentially. 

This heavily controlled noise growth allows the scheme to evaluate relatively deep homomorphic circuits before necessitating bootstrapping.

## Bootstrapping (TFHE / FHEW Blind Rotation)

*Note: While the original project requirements explicitly omitted bootstrapping, I went beyond the scope of the project for learning purposes to implement it here.*

This project includes a fully functional bootstrapping mechanism for Ring-GSW ciphertexts (`src/bootstrapping/mod.rs`), implementing the **AP14 / TFHE (CGGI16) Blind Rotation** algorithm. Bootstrapping clears the accumulated noise of a ciphertext, enabling theoretically infinite homomorphic evaluations (Fully Homomorphic Encryption).

### Mathematical Foundation of Blind Rotation
The core idea is to evaluate the LWE decryption function homomorphically inside the exponent of a polynomial ring $R_q = \mathbb{Z}_q[X]/(X^N + 1)$. The LWE decryption phase is defined as:
$\phi \approx b - \sum a_i s_i \pmod q$

Instead of computing this directly in integers, we rotate a **Test Polynomial** (which acts as a lookup table/step function to map noisy phase values back to clean message boundaries) by this phase degree.
1. **Initialization:** The accumulator (an RGSW encryption of the test polynomial) is uniformly multiplied/rotated by $X^b$.
2. **Conditional Rotations (CMUX):** For each scalar $a_i$ in the LWE ciphertext, we homomorphically rotate the accumulator by $X^{-a_i}$ *if and only if* the secret key bit $s_i = 1$. Since $s_i$ is provided as an RGSW ciphertext inside the **Bootstrapping Key**, it operates as a homomorphic selector (multiplexer):
   $ACC_{new} = ACC_{current} + RGSW(s_i) \cdot (ACC_{current} \cdot X^{-a_i} - ACC_{current})$
3. **Sample Extraction:** After processing all $a_i$ components, the accumulator polynomial has been shifted by exactly the phase $\phi$. The $0$-th degree coefficient of the resulting polynomial naturally falls into the clean region of the lookup table. By extracting the corresponding coefficients from the RGSW matrix, we extract a fresh, low-noise LWE ciphertext.

*Note: The blind rotation inner loop utilizes data-level parallelism via the `rayon` crate to significantly accelerate the heavy polynomial matrix operations.*

You can run the end-to-end integration test demonstrating LWE encryption, modulus switching, blind rotation, and key-switched decryption via:
```bash
cargo test test_integration_lwe_bootstrap -- --nocapture
```

## References

- Craig Gentry, Amit Sahai, and Brent Waters. [Homomorphic encryption from learning with errors: Conceptually-simpler, asymptotically-faster, attribute-based](https://eprint.iacr.org/2013/340.pdf). In Ran Canetti and Juan A. Garay, editors, Advances in Cryptology – CRYPTO 2013, Part I, volume 8042 of Lecture Notes in Computer Science, pages 75–92, Santa Barbara, CA, USA, August 18–22, 2013. Springer, Heidelberg, Germany.
- Jacob Alperin-Sheriff and Chris Peikert. [Faster Bootstrapping with Polynomial Error](https://eprint.iacr.org/2014/233). CRYPTO 2014. (AP14)
- Ilaria Chillotti, Nicolas Gama, Mariya Georgieva, and Malika Izabachène. [Faster Fully Homomorphic Encryption: Bootstrapping in less than 0.1 Seconds](https://eprint.iacr.org/2016/870). ASIACRYPT 2016. (TFHE / CGGI16)