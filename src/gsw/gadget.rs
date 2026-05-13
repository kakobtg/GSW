use nalgebra::DMatrix;
use crate::params::gsw::*;

pub fn gadget_matrix() -> DMatrix<i64> {
    let mut g = DMatrix::zeros(N_LWE + 1, N_COLS);
    for i in 0..=N_LWE {
        for j in 0..L {
            g[(i, i * L + j)] = 1 << j;
        }
    }
    g
}

pub fn g_inverse(c: &DMatrix<i64>) -> DMatrix<i64> {
    let mut z = DMatrix::zeros(N_COLS, N_COLS);
    for col in 0..N_COLS {
        for row in 0..=N_LWE {
            let mut val = c[(row, col)].rem_euclid(Q);
            for j in 0..L {
                z[(row * L + j, col)] = val & 1;
                val >>= 1;
            }
        }
    }
    z
}