# crimlab for JavaScript

Criminology and sociolegal statistics: the C++ kernels of the R and Python
packages [crimlab](https://github.com/rootcoder007/crimlab), built to
WebAssembly, with the `crimlab` command. No dependencies; the module sees
no files and no network.

```sh
npx crimlab orrisk --or 2 --rate 0.1 --share 0.3
echo '{"done": [30, 45, 60], "open": [100, 150]}' | npx crimlab losbnd -j
npx crimlab -l          # every function and its arguments
```

```js
import crimlab, { jrun, jfns } from 'crimlab';

crimlab.orrisk({ or: 2, rate: 0.1, share: 0.3 });
// { exposed: 0.1475..., unexposed: 0.0796..., rr: 1.8524..., or: 2 }
jrun('vodtst', { day: { a: 40, b: 60 }, night: { a: 50, b: 50 }, ref: 'b' });
```

Every function of the Python package is here under its R name, with the
same named arguments, checks and results. A data frame is an object of
columns, a named vector an object of numbers, a missing value `null`, and
a number JSON cannot hold the string `"NaN"`, `"Inf"` or `"-Inf"`. An
object given for a per-group argument names the groups. A refused call
throws an `Error` whose message names the function and what is wrong.

Node 18 or later, or any browser with WebAssembly (the module is fetched
next to `lib/`). Licence: Apache-2.0 WITH LLVM-exception.
