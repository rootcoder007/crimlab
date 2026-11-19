#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""The unity build of the R package (src/Makevars).

R compiles src/u_g*.cpp (the R glue) and src/u_k*.cpp (the kernels and
their self-tests); each includes a group of the sources, so the headers
are read once per group and not once per file. Every other build (the
command line, Python, Go, Rust, JavaScript) takes each source by itself
and passes the unity files by: they include a header of R.

    python3 dev/unity.py check     every source is in exactly one group
    python3 dev/unity.py make [N]  group the sources again, N to a group

`make` starts from groups of N files (22 by default) and moves a file
that does not compile with its group (two files with a helper of one
name in an unnamed namespace) to a later group; it needs g++ and R.
After it, src/Makevars lists the object files it prints.
"""
import concurrent.futures as cf
import os
import re
import subprocess
import sys

top = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = os.path.join(top, "src")
pat = re.compile(r'#\s*include\s*[<"](R\.h|Rinternals\.h|Rmath\.h|R_ext/|rchk\.h)')
head = """// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): %s
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

"""


def sources():
    return sorted(f for f in os.listdir(src) if f.endswith(".cpp") and not f.startswith("u_"))


def units():
    return sorted(f for f in os.listdir(src) if re.fullmatch(r"u_[gk]\d+\.cpp", f))


def check():
    seen = {}
    for u in units():
        for m in re.finditer(r'^#include "([^"]+)"', open(os.path.join(src, u)).read(), re.M):
            seen.setdefault(m.group(1), []).append(u)
    bad = [f for f in sources() if len(seen.get(f, [])) != 1] + [f for f in seen if f not in sources()]
    mk = open(os.path.join(src, "Makevars")).read()
    bad += [u for u in units() if u[:-4] + ".o" not in mk]
    if bad:
        sys.exit("unity.py: not in exactly one group, or not in src/Makevars: " + " ".join(sorted(set(bad))))
    print("unity.py: %d sources in %d groups" % (len(sources()), len(units())))


def fails(group, tag, tmp, rinc):
    p = os.path.join(tmp, "try_%s.cpp" % tag)
    with open(p, "w") as f:
        f.write("".join('#include "%s"\n' % x for x in group))
    r = subprocess.run(["g++", "-std=gnu++17", "-fsyntax-only", "-DCLAB_R", "-I" + src, "-I" + rinc, "-w", p],
                       capture_output=True, text=True)
    if r.returncode == 0:
        return None
    for m in re.finditer(r"/([A-Za-z0-9_]+\.cpp):\d+:\d+: error", r.stderr):
        if m.group(1) in group:
            return m.group(1)
    return group[-1]


def pack(files, tag, per, tmp, rinc):
    groups, rest = [], list(files)
    while rest:
        g, rest = rest[:per], rest[per:]
        while True:
            bad = fails(g, "%s%d" % (tag, len(groups)), tmp, rinc)
            if bad is None:
                break
            g.remove(bad)
            rest.append(bad)
            if not g:
                sys.exit("unity.py: does not compile by itself: " + bad)
        groups.append(g)
    return groups


def make(per):
    import tempfile
    rinc = subprocess.check_output(["R", "--no-echo", "-e", "cat(R.home('include'))"], text=True).strip()
    files = sources()
    glue = [f for f in files if pat.search(open(os.path.join(src, f), errors="replace").read())]
    kern = [f for f in files if f not in glue]
    for u in units():
        os.remove(os.path.join(src, u))
    with tempfile.TemporaryDirectory() as tmp, cf.ThreadPoolExecutor(2) as ex:
        a = ex.submit(pack, glue, "g", per, tmp, rinc)
        b = ex.submit(pack, kern, "k", per, tmp, rinc)
        res = (("g", a.result(), "the R glue"), ("k", b.result(), "the kernels"))
    objs = []
    for tag, groups, what in res:
        for i, g in enumerate(groups, 1):
            nm = "u_%s%02d" % (tag, i)
            objs.append(nm + ".o")
            with open(os.path.join(src, nm + ".cpp"), "w") as f:
                f.write(head % what + "".join('#include "%s"\n' % x for x in g))
    print("OBJECTS = " + " ".join(objs))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "make":
        make(int(sys.argv[2]) if len(sys.argv) > 2 else 22)
    elif len(sys.argv) > 1 and sys.argv[1] == "check":
        check()
    else:
        sys.exit(__doc__)
