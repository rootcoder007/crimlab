#!/usr/bin/env node
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// npx crimlab: the command-line program (cli/main.cpp) built to
// WebAssembly. The module sees no files: -f FILE is read here and given
// to it as standard input.
import { readFileSync } from 'node:fs';
import { load } from '../lib/wasm.js';

const argv = process.argv.slice(2);
let input = null;
const i = argv.indexOf('-f');
if (i > 0 && i + 1 < argv.length && argv[i + 1] !== '-') {
  try {
    input = readFileSync(argv[i + 1]);
  } catch {
    process.stderr.write(`crimlab: cannot open ${argv[i + 1]}.\n`);
    process.exit(2);
  }
  argv[i + 1] = '-';
}
const m = await load(() => input ?? readFileSync(0));
const r = m.cmain(argv);
process.stdout.write(r.out);
process.stderr.write(r.err);
process.exitCode = r.code;
