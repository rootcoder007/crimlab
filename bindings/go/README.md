# crimlab for Go <img src="https://raw.githubusercontent.com/rootcoder007/crimlab/main/man/figures/logo.png" align="right" height="139" alt="crimlab hex logo" />

[![Website](https://img.shields.io/badge/website-crimlab.dev-b8432a)](https://crimlab.dev/)
[![check](https://github.com/rootcoder007/crimlab/actions/workflows/check.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/check.yml)
[![coverage](https://github.com/rootcoder007/crimlab/actions/workflows/coverage.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/coverage.yml)
[![srr](https://github.com/rootcoder007/crimlab/actions/workflows/srr.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/srr.yml)
[![pkgdown](https://github.com/rootcoder007/crimlab/actions/workflows/pkgdown.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/pkgdown.yml)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/16873/badge)](https://www.bestpractices.dev/projects/16873)
[![Release](https://img.shields.io/github/v/release/rootcoder007/crimlab?label=release&color=2ea44f)](https://github.com/rootcoder007/crimlab/releases/latest)
[![PyPI](https://img.shields.io/pypi/v/crimlab?label=PyPI&color=3775a9)](https://pypi.org/project/crimlab/)
[![npm](https://img.shields.io/npm/v/crimlab?label=npm&color=cb3837)](https://www.npmjs.com/package/crimlab)
[![crates.io](https://img.shields.io/crates/v/crimlab?label=crates.io&color=f0b400)](https://docs.rs/crimlab/latest/crimlab/)
[![Maven Central](https://img.shields.io/maven-central/v/io.github.rootcoder007/crimlab?label=Maven%20Central&color=6f42c1)](https://central.sonatype.com/artifact/io.github.rootcoder007/crimlab)
[![Go Reference](https://pkg.go.dev/badge/github.com/rootcoder007/crimlab/bindings/go.svg)](https://pkg.go.dev/github.com/rootcoder007/crimlab/bindings/go)
[![COPR](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fcopr.fedorainfracloud.org%2Fapi_3%2Fpackage%2F%3Fownername%3Drootcoder007%26projectname%3Dcrimlab%26packagename%3Dcrimlab%26with_latest_build%3Dtrue&query=%24.builds.latest.state&label=COPR%20build&color=2ea44f)](https://copr.fedorainfracloud.org/coprs/rootcoder007/crimlab/)
[![OBS](https://build.opensuse.org/projects/home:rootcoder007/packages/crimlab/badge.svg?type=default)](https://build.opensuse.org/package/show/home:rootcoder007/crimlab)
[![r-universe](https://rootcoder007.r-universe.dev/badges/crimlab)](https://rootcoder007.r-universe.dev/crimlab)

Documentation, guides and every function: <https://crimlab.dev/>.

Statistical Toolkit for Criminology and Sociolegal Research: the C++ kernels of the R and Python
packages [crimlab](https://github.com/rootcoder007/crimlab), compiled into
the package with cgo from the copy in `csrc/` (a C++17 compiler is all it
needs). No other dependencies.

Every function of the R package is here with its first letter capitalised
(`orrisk` is `Orrisk`); `Jrun` calls one by its R name. A function takes
its named arguments as a JSON object and returns its result as JSON, with
the names, checks and results of the R and Python packages: a data frame
is an object of columns, a named vector an object of numbers, a missing
value `null`, and a number JSON cannot hold `"NaN"`, `"Inf"` or `"-Inf"`.

```go
package main

import (
	"fmt"

	crimlab "github.com/rootcoder007/crimlab/bindings/go"
)

func main() {
	r, err := crimlab.Orrisk(`{"or": 2, "rate": 0.1, "share": 0.3}`)
	if err != nil {
		panic(err)
	}
	fmt.Println(r)
	// the same by name
	r, err = crimlab.Jrun("losbnd", `{"done": [30, 45, 60], "open": [100, 150]}`)
	fmt.Println(r, err)
}
```

A refused call returns an `*Error`, whose message names the function and
what is wrong, as in R and Python. Licence: Apache-2.0 WITH
LLVM-exception.

## Citation

> Ruhela, V. S. (2026). *crimlab: Statistical Toolkit for Criminology and
> Sociolegal Research.* <https://github.com/rootcoder007/crimlab>

`CITATION.cff` in the repository holds the same reference for reference
managers; in R, `citation("crimlab")` prints it.
