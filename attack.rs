use nalgebra::DMatrix;
use crate::params::gsw::*;

pub fn export_lattice_basis_for_sage(a: &DMatrix<i64>) {
    println!("\n# === SAGE MATH LATTICE ATTACK SCRIPT ===");
    println!("# Copy and paste everything below this line into SageMath (or an online Sage cell) to run the attack:");
    println!("m = {}", M);
    println!("n = {}", N_LWE);
    println!("q = {}", Q);
    println!("basis = [");
    for i in 0..M {
        let mut row = vec![0; M + N_LWE + 1];
        row[i] = Q;
        println!("    {:?},", row);
    }
    for i in 0..N_LWE {
        let mut row = vec![0; M + N_LWE + 1];
        for j in 0..M { row[j] = a[(i, j)]; }
        row[M + i] = 1;
        println!("    {:?},", row);
    }
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
    println!("for row in L_red:\n    if row[-1] == 1 or row[-1] == -1:\n        secret = row[m:m+n]");
    println!("        if row[-1] == -1: secret = [-x for x in secret]");
    println!("        print(\"\\n[+] LWE Broken! Recovered Secret s:\", secret)\nprint(\"\\n[+] FINISHED\")");
}