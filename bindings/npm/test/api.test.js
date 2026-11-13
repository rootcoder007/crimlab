// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { execFileSync, spawnSync } from 'node:child_process';
import crimlab, { jrun, jfns } from '../lib/index.js';

const bin = new URL('../bin/crimlab.js', import.meta.url).pathname;

test('every case of cli/test.tsv, as the Python package gives it', () => {
  const rows = readFileSync(new URL('../../../cli/test.tsv', import.meta.url), 'utf8').trim().split('\n');
  for (const r of rows) {
    const [f, a, want] = r.split('\t');
    const got = jrun(f, a);
    const w = JSON.parse(want);
    const near = (x, y) =>
      typeof x === 'number' && typeof y === 'number'
        ? Math.abs(x - y) <= 1e-12 * Math.abs(y)
        : Array.isArray(x)
          ? x.length === y.length && x.every((v, k) => near(v, y[k]))
          : x && typeof x === 'object'
            ? Object.keys(x).join() === Object.keys(y).join() && Object.keys(x).every((k) => near(x[k], y[k]))
            : x === y;
    assert.ok(near(got, w), `${f} ${a}`);
  }
});

test('methods, errors and the list', () => {
  assert.deepEqual(crimlab.losbnd({ done: [30, 45, 60, 90, 120], open: [100, 150, 200] }),
    { mean: 69, age: 150, lower: 99.375, bias: 30.375 });
  assert.equal(jrun('selftest').nfail, 0);
  assert.deepEqual(jfns().kinidx, ['g1', 'g2', 'freq', 'rel?']);
  assert.throws(() => crimlab.orrisk({ or: 2, rate: 1, share: 0.3 }),
    { message: 'orrisk: `rate` must be one number strictly between 0 and 1.' });
  assert.throws(() => jrun('nope'), /no function nope/);
});

test('npx crimlab', () => {
  const o = execFileSync(process.execPath, [bin, 'losbnd', '-j'], { input: '{"done": [30, 45], "open": [100]}' });
  assert.equal(String(o), '{"mean":37.5,"age":100,"lower":58.333333333333336,"bias":20.833333333333332}\n');
  const t = execFileSync(process.execPath, [bin, 'orrisk', '--or', '2', '--rate', '0.1', '--share', '0.3']);
  assert.match(String(t), /^exposed {5}0\.1475205\n/);
  const e = spawnSync(process.execPath, [bin, 'losbnd', '--done', '1'], { input: '' });
  assert.equal(e.status, 1);
  assert.equal(String(e.stderr), 'losbnd: `open` is missing.\n');
  assert.equal(String(execFileSync(process.execPath, [bin, '-V'])).trim().split(' ')[0], 'crimlab');
});
