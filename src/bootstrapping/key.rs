use rand::Rng;
use crate::rgsw::{encrypt, PolyMatrix, RgswBitCiphertext};

#[derive(Clone, Debug)]
pub struct BootstrappingKey {
    pub bk: Vec<RgswBitCiphertext>,
}

pub fn generate_bootstrapping_key(rgsw_pk: &PolyMatrix, lwe_sk: &[i64], rng: &mut impl Rng) -> BootstrappingKey {
    let mut bk = Vec::new();
    for &sk_bit in lwe_sk {
        bk.push(encrypt(rgsw_pk, sk_bit, rng).unwrap());
    }
    BootstrappingKey { bk }
}