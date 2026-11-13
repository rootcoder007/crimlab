# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Fedora COPR, build method "make_srpm": the source RPM from this
# checkout, its version from DESCRIPTION.
#   make -f .copr/Makefile srpm outdir=DIR
outdir ?= $(CURDIR)
top := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
ver := $(shell sed -n 's/^Version: *//p' $(top)/DESCRIPTION)
tmp := $(outdir)/.crimlab-srpm

srpm:
	command -v git >/dev/null || dnf -y install git
	rm -rf $(tmp) && mkdir -p $(tmp)
	git -C $(top) config --global --add safe.directory $(top)
	git -C $(top) archive --prefix=crimlab-$(ver)/ -o $(tmp)/crimlab-$(ver).tar.gz HEAD
	sed 's/^Version:.*/Version:        $(ver)/' $(top)/packaging/rpm/crimlab.spec > $(tmp)/crimlab.spec
	rpmbuild -bs --define "_sourcedir $(tmp)" --define "_srcrpmdir $(outdir)" $(tmp)/crimlab.spec
	rm -rf $(tmp)

.PHONY: srpm
