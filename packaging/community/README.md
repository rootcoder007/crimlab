# Community recipes

The files each outside repository wants when crimlab is submitted there.
They build the published release tarball and carry its checksum.

| Directory | Goes to |
|---|---|
| `alpine/testing/crimlab` | gitlab.alpinelinux.org alpine/aports |
| `conan/recipes/crimlab` | conan-io/conan-center-index |
| `conda-forge/recipes/crimlab` | conda-forge/staged-recipes |
| `t2/package/scientific/crimlab` | T2 SDE (svn.exactcode.de t2/trunk, mirrored to github.com/rxrbln/t2sde) |
| `winget/manifests` | microsoft/winget-pkgs |

`.github/workflows/pkgs.yml` builds and tests every one of them on its own
system when they change and when a release is published, and builds this
tree on glibc 2.17.

At a release the version, the URL and the checksums here move to the new
tarball once it is published: the recipes name a tarball that exists.
The T2 checksum is T2's own, the SHA-224 of the decompressed tarball
(`gzip -dc crimlab-X.Y.Z.tar.gz | sha224sum`), and `crimlab.cache` is the
file a T2 build of the package writes (`t2/ci.sh` prints it).
