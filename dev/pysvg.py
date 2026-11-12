# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Writes py/tipplt.svg, the picture in py/README.md, from the README's own
example: the tipplt() call in its Plots section, read from the README and
run with the installed package. From the repository root:
    python dev/pysvg.py
"""

import ast
import re
from pathlib import Path

import crimlab

readme = Path("py/README.md").read_text(encoding="utf-8")
inner = re.search(r">>> svg = tipplt\((.*?)\)\.svg\(\)", readme).group(1)
same, diff = ast.literal_eval(f"({inner})")  # literals only, never code
svg = crimlab.tipplt(same, diff).svg()
Path("py/tipplt.svg").write_text(svg, encoding="utf-8")
print(f"py/tipplt.svg: {len(svg)} bytes from tipplt({inner})")
