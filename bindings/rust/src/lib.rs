// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//! crimlab: criminology and sociolegal statistics.
//!
//! The C++ kernels of the R and Python packages crimlab, compiled into
//! the crate by its build script. Every function takes its named
//! arguments as a JSON object and returns its result as JSON, with the
//! names, checks and results of the R and Python packages: a data frame
//! is an object of columns, a named vector an object of numbers, a
//! missing value `null`, and a number JSON cannot hold `"NaN"`, `"Inf"`
//! or `"-Inf"`.
//!
//! ```
//! let r = crimlab::losbnd(r#"{"done": [30, 45, 60, 90, 120], "open": [100, 150, 200]}"#)?;
//! assert_eq!(r, r#"{"mean":69,"age":150,"lower":99.375,"bias":30.375}"#);
//! # Ok::<(), crimlab::Error>(())
//! ```
#![deny(missing_docs)]

use std::ffi::{CStr, CString};
use std::fmt;
use std::os::raw::{c_char, c_int};

extern "C" {
    fn clab_jrun(f: *const c_char, a: *const c_char, out: *mut c_char, cap: c_int, need: *mut c_int) -> c_int;
    fn clab_jfns() -> *const c_char;
}

/// A refused call: the message names the function and what is wrong.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

/// A JSON result, or the message of a refused call.
pub type Result<T> = std::result::Result<T, Error>;

/// Runs the function named `name` (an R function name) on `args`, a JSON
/// object of named arguments.
pub fn jrun(name: &str, args: &str) -> Result<String> {
    let f = CString::new(name).map_err(|_| Error("crimlab: a NUL in the function name.".into()))?;
    let a = CString::new(args).map_err(|_| Error("crimlab: a NUL in the arguments.".into()))?;
    let mut buf = vec![0u8; 1 << 16];
    loop {
        let mut need: c_int = 0;
        let cap = c_int::try_from(buf.len()).unwrap_or(c_int::MAX);
        // SAFETY: f and a are NUL-terminated; buf holds cap bytes
        let rc = unsafe { clab_jrun(f.as_ptr(), a.as_ptr(), buf.as_mut_ptr().cast(), cap, &mut need) };
        if rc == 3 && need as usize > buf.len() {
            buf.resize(need as usize, 0);
            continue;
        }
        let s = CStr::from_bytes_until_nul(&buf).unwrap_or_default().to_string_lossy().into_owned();
        return if rc == 0 { Ok(s) } else { Err(Error(s)) };
    }
}

/// The functions and their arguments, as a JSON object from name to
/// argument names (`?` marks an optional one).
pub fn jfns() -> &'static str {
    // SAFETY: a NUL-terminated string in static storage
    unsafe { CStr::from_ptr(clab_jfns()) }.to_str().unwrap_or("{}")
}

include!(concat!(env!("OUT_DIR"), "/fns.rs"));
