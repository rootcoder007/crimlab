# Security policy

crimlab is a statistical library: it reads the data a caller passes, and
opens no network connection except in the documented readers of open
data
([`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md),
[`odls()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)).
It implements no cryptography. It has not had an external security
review.

Please report a vulnerability privately through GitHub’s “Report a
vulnerability” button on this repository (Security tab,
<https://github.com/rootcoder007/crimlab/security/advisories/new>), not
in a public issue. A reply comes within 7 days, and a fix and an
advisory within 90 days. Fixed vulnerabilities are listed at
<https://github.com/rootcoder007/crimlab/security/advisories>.

A crash, a hang or memory that grows without bound on input a caller
could reasonably pass counts as a vulnerability here, in the R package,
the command-line program and every binding. A wrong number is a bug:
report it in the issue tracker.
