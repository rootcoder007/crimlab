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

## Release signing key

Each release’s `SHA256SUMS` is signed with this SSH key (namespace
`file`); the signature is the release file `SHA256SUMS.sig`.

    ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIDCgdm/W6LtsPET14oTqb1lOFga52VG4HSJTtQHfLTFc

Fingerprint `SHA256:o8zFxDLZk0Q4Cm00HWgHXqfD/rlI6TVLAI00OaAtrAk`. To
check a release, put the key in an allowed-signers file and verify the
sums, then check the archive against them:

``` sh
echo "crimlab ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIDCgdm/W6LtsPET14oTqb1lOFga52VG4HSJTtQHfLTFc" > allowed_signers
ssh-keygen -Y verify -f allowed_signers -I crimlab -n file \
  -s SHA256SUMS.sig < SHA256SUMS
sha256sum -c --ignore-missing SHA256SUMS
```

`install.sh` and `install.ps1` make the same check when the environment
variable `CRIMLAB_SIGNER` holds that allowed-signers line
(`crimlab ssh-ed25519 ...`); without it they check the archive against
`SHA256SUMS` only.
