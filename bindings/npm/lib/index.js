// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// crimlab: criminology and sociolegal statistics. Every function of the R
// and Python packages crimlab, by its R name, from the C++ kernels built
// to WebAssembly. A function takes its named arguments as an object and
// returns its result as JSON.parse gives it: a data frame is an object of
// columns, a named vector an object of numbers, a missing value null, and
// a number JSON cannot hold the string "NaN", "Inf" or "-Inf".
import { load } from './wasm.js';

const m = await load();

/** Runs the function `name` on `args` (an object, or JSON text). */
export function jrun(name, args = {}) {
  const [rc, s] = m.jrun(String(name), typeof args === 'string' ? args : JSON.stringify(args));
  if (rc !== 0) throw new Error(s);
  return JSON.parse(s);
}

/** The functions and their argument names ("?" marks an optional one). */
export function jfns() {
  return JSON.parse(m.jfns());
}

/** Every function as a method: crimlab.orrisk({ or: 2, rate: 0.1, share: 0.3 }). */
const crimlab = Object.freeze(
  Object.fromEntries(Object.keys(jfns()).map((n) => [n, (args) => jrun(n, args)])),
);
export default crimlab;
