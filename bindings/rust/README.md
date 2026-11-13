# crimlab for Rust

Criminology and sociolegal statistics: the C++ kernels of the R and Python
packages [crimlab](https://github.com/rootcoder007/crimlab), compiled into
the crate by its build script with the `cc` crate (a C++17 compiler is all
it needs). No other dependencies.

Every function of the Python package is here under its R name. It takes
its named arguments as a JSON object and returns its result as JSON, with
the names, checks and results of the R and Python packages: a data frame
is an object of columns, a named vector an object of numbers, a missing
value `null`, and a number JSON cannot hold `"NaN"`, `"Inf"` or `"-Inf"`.
A JSON object given for a per-group argument names the groups.

```rust
fn main() -> crimlab::Result<()> {
    let r = crimlab::orrisk(r#"{"or": 2, "rate": 0.1, "share": 0.3}"#)?;
    println!("{r}");
    // the same by name; crimlab::jfns() lists every function
    let r = crimlab::jrun("losbnd", r#"{"done": [30, 45, 60], "open": [100, 150]}"#)?;
    println!("{r}");
    Ok(())
}
```

A refused call returns `crimlab::Error`, whose message names the function
and what is wrong, as in R and Python. Licence: Apache-2.0 WITH
LLVM-exception.
