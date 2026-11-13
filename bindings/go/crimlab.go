// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

// Package crimlab is the Go binding of crimlab, criminology and
// sociolegal statistics: the C++ kernels of the R and Python packages
// crimlab, compiled into the package with cgo from the copy in csrc/
// (tools/vendor.sh makes it); a C++17 compiler is all it needs.
//
// Every function takes its named arguments as a JSON object and returns
// its result as JSON, with the names, checks and results of the R and
// Python packages: a data frame is an object of columns, a named vector
// an object of numbers, a missing value null, and a number JSON cannot
// hold "NaN", "Inf" or "-Inf". The functions are those of the R package,
// with the first letter capitalised (orrisk is Orrisk); Jrun calls one
// by its R name.
package crimlab

/*
#cgo CXXFLAGS: -std=c++17 -O2 -I${SRCDIR}/csrc/src -I${SRCDIR}/csrc/cli
#include <stdlib.h>
int clab_jrun(const char *fn, const char *in, char *out, int cap, int *need);
const char *clab_jfns(void);
*/
import "C"

import (
	"strings"
	"unsafe"
)

// Error is a refused call: Msg names the function and what is wrong.
type Error struct{ Msg string }

func (e *Error) Error() string { return e.Msg }

// Jrun runs the function named name (an R function name, such as
// "orrisk") on args, a JSON object of named arguments ("{}" or "" for
// none), and returns its JSON result.
func Jrun(name, args string) (string, error) {
	if strings.IndexByte(name, 0) >= 0 || strings.IndexByte(args, 0) >= 0 {
		return "", &Error{"crimlab: a NUL in the function name or the arguments."}
	}
	if args == "" {
		args = "{}"
	}
	f := C.CString(name)
	defer C.free(unsafe.Pointer(f))
	a := C.CString(args)
	defer C.free(unsafe.Pointer(a))
	n := 1 << 16
	for {
		buf := (*C.char)(C.malloc(C.size_t(n)))
		var need C.int
		rc := C.clab_jrun(f, a, buf, C.int(n), &need)
		if rc == 3 && int(need) > n { // CLAB_ESIZE: again, with room
			C.free(unsafe.Pointer(buf))
			n = int(need)
			continue
		}
		s := C.GoString(buf)
		C.free(unsafe.Pointer(buf))
		if rc != 0 {
			return "", &Error{s}
		}
		return s, nil
	}
}

// Jfns returns the functions and their arguments, as a JSON object from
// name to argument names ("?" marks an optional one).
func Jfns() string { return C.GoString(C.clab_jfns()) }
