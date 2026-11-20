#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""The release notes of one version: its section of NEWS.md with each
paragraph and each list item on one line. NEWS.md is wrapped at 80 columns
for R; GitHub shows a release body as written, so a wrapped one reads as
broken lines.
  python3 dev/notes.py NEWS.md 0.0.4 > notes.md
Exits 1 when NEWS.md has no section "# crimlab <version>"."""
import re
import sys


def section(text, version):
    out, on = [], False
    for line in text.replace("\r\n", "\n").split("\n"):
        if line == f"# crimlab {version}":
            on = True
            continue
        if on and line.startswith("# "):
            break
        if on:
            out.append(line)
    return out if on else None


def unwrap(lines):
    out, fence = [], False
    for line in lines:
        if line.strip().startswith("```"):
            fence = not fence
            out.append(line)
            continue
        starts = (
            fence
            or not line.strip()
            or not out
            or not out[-1].strip()
            or out[-1].strip().startswith("```")
            or re.match(r"\s*([-*+] |\d+[.)] |#|\||>)", line)
            or out[-1].lstrip().startswith(("#", "|"))
        )
        if starts:
            out.append(line)
        else:
            out[-1] = out[-1].rstrip() + " " + line.strip()
    return "\n".join(out).strip("\n") + "\n"


def main():
    news, version = sys.argv[1], sys.argv[2]
    with open(news, encoding="utf-8") as f:
        lines = section(f.read(), version)
    if lines is None or not any(x.strip() for x in lines):
        sys.exit(f"notes.py: {news} has no section for version {version}")
    sys.stdout.write(unwrap(lines))


if __name__ == "__main__":
    main()
