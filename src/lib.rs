pub mod attack;
pub mod bootstrapping;
pub mod gsw;
pub mod params;
pub mod rgsw;

#[derive(Debug)]
pub enum FheError {
    InvalidMessage { value: i64, expected: &'static str },
    DimensionMismatch { got: usize, expected: usize },
    ValueOutOfRange { value: i64, max_abs: i64 },
    VectorLengthMismatch { len1: usize, len2: usize },
}