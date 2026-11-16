// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The host side of crimlab.wasm (build.sh): the WASI preview 1 functions
// it imports. It sees no files, no environment and no network: standard
// input is the bytes the caller gives (or, for the command, this
// process's standard input, read when the module first asks), and what it
// writes to standard output and error is collected. Every other call
// answers EBADF or ENOSYS.

const BADF = 8;
const NOENT = 44;
const NOSYS = 52;
const SPIPE = 70;
const enc = new TextEncoder();
const dec = new TextDecoder();

class Exit {
  constructor(code) {
    this.code = code;
  }
}

async function bytes() {
  const u = new URL('../dist/crimlab.wasm', import.meta.url);
  if (u.protocol === 'file:') {
    const { readFile } = await import('node:fs/promises');
    return readFile(u);
  }
  const r = await fetch(u);
  if (!r.ok) throw new Error(`crimlab: cannot fetch ${u}`);
  return new Uint8Array(await r.arrayBuffer());
}

/** Instantiates the module. stdin() returns the bytes of standard input. */
export async function load(stdin = () => new Uint8Array(0)) {
  let mem = null;
  let inp = null;
  let inpos = 0;
  const out = [[], []];
  const dv = () => new DataView(mem.buffer);
  const u8 = () => new Uint8Array(mem.buffer);
  const iov = (p, n, f) => {
    let t = 0;
    for (let i = 0; i < n; i++) {
      const b = dv().getUint32(p + 8 * i, true);
      const l = dv().getUint32(p + 8 * i + 4, true);
      const k = f(b, l);
      t += k;
      if (k < l) break;
    }
    return t;
  };
  const wasi = {
    // an empty environment: no variables, no bytes (the C++ library asks when it sets up its locale)
    environ_sizes_get(np, sp) {
      dv().setUint32(np, 0, true);
      dv().setUint32(sp, 0, true);
      return 0;
    },
    environ_get: () => 0,
    clock_time_get(id, prec, p) {
      dv().setBigUint64(p, BigInt(Math.round(performance.now() * 1e6)), true);
      return 0;
    },
    fd_close: () => 0,
    fd_fdstat_get(fd, p) {
      if (fd > 2) return BADF;
      u8().fill(0, p, p + 24);
      dv().setUint8(p, 2); // a character device
      return 0;
    },
    fd_fdstat_set_flags: () => NOSYS,
    fd_prestat_get: () => BADF,
    fd_prestat_dir_name: () => BADF,
    path_open: () => NOENT,
    fd_seek: () => SPIPE,
    fd_read(fd, p, n, np) {
      if (fd !== 0) return BADF;
      if (inp === null) inp = stdin();
      dv().setUint32(np, iov(p, n, (b, l) => {
        const k = Math.min(l, inp.length - inpos);
        u8().set(inp.subarray(inpos, inpos + k), b);
        inpos += k;
        return k;
      }), true);
      return 0;
    },
    fd_write(fd, p, n, np) {
      if (fd !== 1 && fd !== 2) return BADF;
      dv().setUint32(np, iov(p, n, (b, l) => {
        out[fd - 1].push(u8().slice(b, b + l));
        return l;
      }), true);
      return 0;
    },
    proc_exit(code) {
      throw new Exit(code);
    },
  };
  // There is no file system here. A function that would read a folder or a
  // file (the power counters of raplok(), the time-zone files of dtctloc())
  // imports calls this object does not list; each answers "not supported",
  // so the function refuses or falls back as it does where the file is
  // absent, and a new import cannot stop the module from loading.
  const host = new Proxy(wasi, { get: (t, k) => (k in t ? t[k] : () => NOSYS) });
  const { instance } = await WebAssembly.instantiate(await bytes(), { wasi_snapshot_preview1: host });
  const x = instance.exports;
  mem = x.memory;
  x._initialize();

  const cstr = (s) => {
    const b = enc.encode(s);
    const p = x.malloc(b.length + 1);
    u8().set(b, p);
    u8()[p + b.length] = 0;
    return p;
  };
  const read = (p) => {
    const m = u8();
    let e = p;
    while (m[e] !== 0) e++;
    return dec.decode(m.subarray(p, e));
  };
  const drain = (i) => {
    const s = dec.decode(Buffer.concat(out[i]));
    out[i] = [];
    return s;
  };

  // A failed allocation cannot be unwound in a build without exceptions:
  // the module's new-handler records it (clab_oomed()) and stops it. Such a
  // stop is reported in words, as the other bindings report CLAB_ENOMEM,
  // and the module, whose heap is then unusable, refuses further calls.
  let dead = false;
  const oom = (e) => e instanceof WebAssembly.RuntimeError && x.clab_oomed && x.clab_oomed() !== 0;
  const alive = () => {
    if (dead) throw new Error('crimlab: the module ran out of memory; load it again');
  };

  return {
    /** clab_jrun(): [status, JSON result or message]. */
    jrun(name, args) {
      alive();
      let f, a, need;
      try {
        f = cstr(name);
        a = cstr(args);
        need = x.malloc(4);
      } catch (e) {
        if (!oom(e)) throw e;
        dead = true;
        return [4, `${name}: not enough memory for this size.`];
      }
      const a = cstr(args);
      const need = x.malloc(4);
      let cap = 1 << 16;
      try {
        for (;;) {
          const o = x.malloc(cap);
          const rc = x.clab_jrun(f, a, o, cap, need);
          const n = dv().getInt32(need, true);
          const s = rc === 3 && n > cap ? null : read(o);
          x.free(o);
          if (s === null) {
            cap = n;
            continue;
          }
          return [rc, s];
        }
      } catch (e) {
        if (!oom(e)) throw e;
        dead = true;
        return [4, `${name}: not enough memory for this size.`];
      } finally {
        if (!dead) {
          x.free(f);
          x.free(a);
          x.free(need);
        }
      }
    },
    jfns: () => read(x.clab_jfns()),
    /** The command: its exit status, standard output and standard error. */
    cmain(argv) {
      alive();
      const ps = ['crimlab', ...argv].map(cstr);
      const v = x.malloc(4 * (ps.length + 1));
      ps.forEach((p, i) => dv().setUint32(v + 4 * i, p, true));
      dv().setUint32(v + 4 * ps.length, 0, true);
      let code;
      try {
        code = x.clab_cmain(ps.length, v);
      } catch (e) {
        if (oom(e)) {
          dead = true;
          return { code: 1, out: drain(0), err: drain(1) + 'crimlab: not enough memory for this size.\n' };
        }
        if (!(e instanceof Exit)) throw e;
        code = e.code;
      }
      x.fflush(0);
      return { code, out: drain(0), err: drain(1) };
    },
  };
}
