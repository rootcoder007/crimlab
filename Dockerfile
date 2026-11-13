# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
# ghcr.io/rootcoder007/crimlab: the crimlab command on Wolfi. The build
# stage compiles cli/ and runs its tests; the final stage holds the
# program (linked statically, as in the release archives), its manual
# page and the licence files.
#   docker build -t crimlab .        (or podman build)
#   docker run --rm crimlab crimlab orrisk --or 2 --rate 0.1 --share 0.3
#   docker run --rm -i -v "$PWD:/w" -w /w crimlab crimlab losbnd -f args.json
# The base image is pinned by digest; the tag is in the comment.

# cgr.dev/chainguard/wolfi-base:latest
FROM cgr.dev/chainguard/wolfi-base@sha256:0b179206ab162b12a43bacea1b525a4023e0272005c92d08555c61b1beaedc40 AS build
RUN apk add --no-cache build-base
WORKDIR /src
COPY . .
RUN make -C cli -j"$(nproc)" crimlab CXXFLAGS=-O2 LDFLAGS=-static WERROR= \
 && strip cli/crimlab \
 && make -C cli check

# cgr.dev/chainguard/wolfi-base:latest
FROM cgr.dev/chainguard/wolfi-base@sha256:0b179206ab162b12a43bacea1b525a4023e0272005c92d08555c61b1beaedc40
COPY --from=build /src/cli/crimlab /usr/bin/crimlab
COPY --from=build /src/cli/crimlab.1 /usr/share/man/man1/crimlab.1
COPY LICENSE NOTICE /usr/share/licenses/crimlab/
USER 65532
CMD ["crimlab", "-h"]
LABEL org.opencontainers.image.title="crimlab" \
      org.opencontainers.image.description="The crimlab command: criminology and sociolegal statistics" \
      org.opencontainers.image.licenses="Apache-2.0 WITH LLVM-exception" \
      org.opencontainers.image.source="https://github.com/rootcoder007/crimlab"
