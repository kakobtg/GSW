# GSW Homomorphic Encryption in Rust

A functional implementation of the Gentry-Sahai-Waters (GSW) homomorphic encryption scheme (GSW13) in Rust. This project demonstrates Key Generation, Encryption, Decryption, and Homomorphic operations (Addition and Multiplication) based on the Learning with Errors (LWE) problem. It also includes an implementation of **Ring-GSW (RGSW)** based on Ring-LWE.

## Prerequisites

- [Rust](https://www.rust-lang.org/tools/install) (includes Cargo)

The project relies on the following crates (defined in `Cargo.toml`):
- `nalgebra` for matrix operations.
- `rand` for random number generation.

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

This project also features an RGSW implementation (`src/rgsw.rs`) which transitions from standard LWE to Ring-LWE. 

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

## Running the Project

To compile and run the test script, execute the following command from the root of the project directory:

```bash
cargo run
```

This will output the results of encrypting bits, decrypting them, and evaluating homomorphic addition and multiplication (which act as XOR and AND gates, respectively).

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

## Constraints

As per the project requirements, bootstrapping is intentionally omitted from this implementation.

## References

- Craig Gentry, Amit Sahai, and Brent Waters. [Homomorphic encryption from learning with errors: Conceptually-simpler, asymptotically-faster, attribute-based](https://eprint.iacr.org/2013/340.pdf). In Ran Canetti and Juan A. Garay, editors, Advances in Cryptology – CRYPTO 2013, Part I, volume 8042 of Lecture Notes in Computer Science, pages 75–92, Santa Barbara, CA, USA, August 18–22, 2013. Springer, Heidelberg, Germany.