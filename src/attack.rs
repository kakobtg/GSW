use nalgebra::DMatrix;
use crate::params::gsw::*;
use std::fs::File;
use std::io::Write;

pub fn export_lattice_basis_for_sage(a: &DMatrix<i64>) {
    let filename = "lattice_attack.py";
    let mut file = File::create(filename).expect("Failed to create lattice attack file");

    writeln!(file, "# === SAGE MATH LATTICE ATTACK SCRIPT ===").unwrap();
    writeln!(file, "from sage.all import *").unwrap();
    writeln!(file, "m = {}", M).unwrap();
    writeln!(file, "n = {}", N_LWE).unwrap();
    writeln!(file, "q = {}", Q).unwrap();
    writeln!(file, "basis = [").unwrap();
    for i in 0..M {
        let mut row = vec![0; M + N_LWE + 1];
        row[i] = Q;
        writeln!(file, "    {:?},", row).unwrap();
    }
    for i in 0..N_LWE {
        let mut row = vec![0; M + N_LWE + 1];
        for j in 0..M { row[j] = a[(i, j)]; }
        row[M + i] = 1;
        writeln!(file, "    {:?},", row).unwrap();
    }
    {
        let mut row = vec![0; M + N_LWE + 1];
        for j in 0..M { row[j] = -a[(N_LWE, j)]; }
        row[M + N_LWE] = 1;
        writeln!(file, "    {:?}", row).unwrap(); 
    }
    writeln!(file, "]").unwrap();
    writeln!(file, "L = Matrix(ZZ, basis)").unwrap();
    writeln!(file, "print(\"Running LLL / BKZ lattice reduction...\")").unwrap();
    writeln!(file, "L_red = L.LLL() # Or use L.BKZ(block_size=10) for stronger reduction").unwrap();
    writeln!(file, "for row in L_red:\n    if row[-1] == 1 or row[-1] == -1:\n        secret = row[m:m+n]").unwrap();
    writeln!(file, "        if row[-1] == -1: secret = [-x for x in secret]").unwrap();
    writeln!(file, "        print(\"\\n[+] LWE Broken! Recovered Secret s:\", secret)\nprint(\"\\n[+] FINISHED\")").unwrap();

    println!("    -> Successfully exported lattice basis to `{}`!", filename);
    println!("    -> To execute the attack, run: sage {}", filename);
}