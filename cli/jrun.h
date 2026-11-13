/* SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
 *
 * Every function of the Python package by name, with JSON in and JSON out:
 * the one entry point the command-line program and the Rust, Java and
 * JavaScript bindings call. The functions, their arguments and their
 * checks are those of the Python package (py/ext/clab.cpp), so a binding
 * holds no per-function code.
 */
#ifndef CLAB_JRUN_H
#define CLAB_JRUN_H

#include "crimlab.h"

#ifdef __cplusplus
extern "C" {
#endif

/* out was too small: *need bytes (with the final NUL) are needed */
#define CLAB_ESIZE 3

/* Runs function fn (an R function name, such as "orrisk") on in, a JSON
 * object of named arguments ("{}" or NULL for none). Writes to out, of
 * cap bytes, a NUL-terminated JSON result (CLAB_OK) or a one-line message
 * "fn: what is wrong" (CLAB_EARG); *need is the size either needs. A
 * number the JSON cannot hold is written as the string "NaN", "Inf" or
 * "-Inf"; a missing value (R's NA) as null. Thread-safe. */
int clab_jrun(const char *fn, const char *in, char *out, int cap, int *need);

/* The functions and their arguments, as a JSON object
 * {"orrisk": ["or", "rate", "share"], ...}. Static storage. */
const char *clab_jfns(void);

#ifdef __cplusplus
}
#endif

#endif
