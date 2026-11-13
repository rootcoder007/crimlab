// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//! Compiles the crimlab kernels into the crate: every .cpp of src/ that
//! does not include R's headers (the rule py/setup.py follows) and
//! cli/jrun.cpp, from the repository when the crate sits in it, else from
//! vendor/ (tools/vendor.sh). Writes one Rust function per entry of the
//! table in jrun.cpp, so the functions here are those of clab_jrun().

use std::fs;
use std::path::{Path, PathBuf};

fn dirs() -> (PathBuf, PathBuf) {
    let me = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let top = me.join("..").join("..");
    if top.join("src").join("crimlab.h").is_file() && top.join("cli").join("jrun.cpp").is_file() {
        return (top.join("src"), top.join("cli"));
    }
    let v = me.join("vendor");
    if !v.join("src").join("crimlab.h").is_file() {
        panic!("crimlab: vendor/ is missing and ../../src is not here (run tools/vendor.sh)");
    }
    (v.join("src"), v.join("cli"))
}

fn rbound(f: &Path) -> bool {
    let s = fs::read_to_string(f).unwrap();
    s.lines().any(|l| {
        let t: String = l.chars().filter(|c| !c.is_whitespace()).collect();
        ["R.h", "Rinternals.h", "Rmath.h", "R_ext/", "rchk.h"]
            .iter()
            .any(|h| t.starts_with(&format!("#include<{h}")) || t.starts_with(&format!("#include\"{h}")))
    })
}

fn main() {
    let (src, cli) = dirs();
    let mut files: Vec<PathBuf> = fs::read_dir(&src)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "cpp") && !rbound(p))
        .collect();
    files.sort();
    files.push(cli.join("jrun.cpp"));
    let mut b = cc::Build::new();
    b.cpp(true).include(&src).include(&cli).warnings(false).files(&files);
    if b.get_compiler().is_like_msvc() {
        b.flag("/std:c++17").flag("/EHsc");
    } else {
        b.flag("-std=c++17");
    }
    b.compile("crimlab");

    // the table of jrun.cpp: {"name", "args", f},
    let t = fs::read_to_string(cli.join("jrun.cpp")).unwrap();
    let a = t.find("const Ent ents[] = {").expect("crimlab: no table in jrun.cpp");
    let mut out = String::from("// Made by build.rs from the table of cli/jrun.cpp.\n");
    for l in t[a..].lines().skip(1).take_while(|l| !l.starts_with("};")) {
        let q: Vec<&str> = l.split('"').collect();
        if q.len() < 5 {
            continue;
        }
        let (name, args) = (q[1], q[3]);
        let doc = if args.is_empty() { "no arguments".to_string() } else { format!("arguments: {args}") };
        out += &format!(
            "/// `{name}()` of the R package ({doc}; `?` marks an optional one).\n\
             pub fn {name}(args: &str) -> Result<String> {{ jrun(\"{name}\", args) }}\n"
        );
    }
    let o = PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("fns.rs");
    fs::write(o, out).unwrap();
    println!("cargo:rerun-if-changed={}", src.display());
    println!("cargo:rerun-if-changed={}", cli.display());
}
