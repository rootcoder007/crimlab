#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Succeeds when the commit has a green run of the workflow copr into the
# scratch project with every target: the Fedora packages were built on
# COPR's own builders, riscv64 and aarch64 among them, before the tag.
#   sh dev/coprgate.sh OWNER/REPO SHA      (gh, logged in or GH_TOKEN)
set -eu
n=$(gh api "repos/$1/actions/workflows/copr.yml/runs?head_sha=$2&status=success&per_page=100" \
  -q '[.workflow_runs[] | select(.display_title == "copr: crimlab-ci every target")] | length')
if [ "$n" -gt 0 ]; then
  echo "coprgate: $n green scratch build(s) of $2"
else
  echo "coprgate: no green run of the workflow copr (crimlab-ci, every target) for $2" >&2
  echo "coprgate: start it on this commit and tag when it is green" >&2
  exit 1
fi
