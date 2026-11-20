# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
# The command-line program crimlab and the library libcrimlab: the
# kernels of ../src (every .cpp there that does not include R's headers,
# the rule py/setup.py follows) and clab_jrun(). GNU make and any C++17
# compiler; no other tool, no download.
#   make [CXX=clang++]     crimlab, libcrimlab.so.0 (the program links
#                          libcrimlab.a, so it needs no library at run time)
#   make check             the tests (test.sh)
#   make install DESTDIR=/tmp/x PREFIX=/usr LIBDIR=/usr/lib64
# Distribution builds pass WERROR= so a new compiler warning does not
# stop them; CI keeps -Werror.
CXXFLAGS ?= -O2 -g
# 32-bit x86: the arithmetic of SSE2, as on every other platform. The x87
# unit computes in 80 bits and rounds a second time when a value is stored,
# so results differ in their last digits with the compiler's options (the
# openSUSE i586 build failed three checks of the table that way). A build
# that names -mfpmath itself keeps its choice: dev/port.sh tests x87 so.
MACH := $(shell $(CXX) -dumpmachine 2>/dev/null)
ifneq ($(filter i386-% i486-% i586-% i686-%,$(MACH))$(findstring -m32,$(CXX) $(CXXFLAGS)),)
ifeq ($(findstring -mfpmath,$(CXX) $(CXXFLAGS)),)
override CXXFLAGS += -msse2 -mfpmath=sse
endif
endif
WERROR ?= -Werror
WARN = -std=c++17 -Wall -Wextra -pedantic $(WERROR)
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include
MANDIR ?= $(PREFIX)/share/man
SOVER = 0
VER := $(shell sed -n 's/^Version: *//p' ../DESCRIPTION)
# src/ver.h carries the same version for the builds made without this file
ifeq ($(shell grep -c '"$(VER)"' ../src/ver.h),0)
$(error src/ver.h does not name version $(VER) of DESCRIPTION)
endif
B = build

R = '\#[[:space:]]*include[[:space:]]*[<"](R\.h|Rinternals\.h|Rmath\.h|R_ext/|rchk\.h)'
KS := $(shell grep -L -E $(R) ../src/*.cpp)
KO := $(patsubst ../src/%.cpp,$(B)/%.o,$(KS)) $(B)/jrun.o

ifeq ($(shell uname -s),Darwin)
SO = libcrimlab.$(SOVER).dylib
SOFLAGS = -dynamiclib -install_name $(LIBDIR)/$(SO)
else
SO = libcrimlab.so.$(SOVER)
SOFLAGS = -shared -Wl,-soname,$(SO)
endif

# The libraries of the loader (dlopen) and of the threads, named after the
# objects so that --as-needed keeps them: glibc before 2.34 has both outside
# libc, and a link without them stops at dlopen and dlsym there.
# By the machine the compiler builds for (MACH), so a cross build links
# what its target has: macOS and Windows keep both in their system library.
ifneq ($(findstring linux,$(MACH)),)
LDLIBS ?= -pthread -ldl
else ifneq ($(findstring mingw,$(MACH))$(findstring windows,$(MACH)),)
# WinHTTP, a part of Windows: the program asks through it whether a newer
# release exists
LDLIBS ?= -lwinhttp
else ifeq ($(findstring darwin,$(MACH)),)
LDLIBS ?= -pthread
endif

all: crimlab libcrimlab.a $(SO)

$(B)/%.o: ../src/%.cpp ../src/*.h | $(B)
	$(CXX) $(WARN) $(CPPFLAGS) $(CXXFLAGS) -fPIC -c $< -o $@

# The JSON layer is one large file, the costliest to compile by far. JOPT
# is added to its flags alone, after the others: packaging/rpm lowers its
# optimisation and its debug information on riscv64, where the builders
# are emulated. The kernels keep the flags they are given.
$(B)/jrun.o: jrun.cpp jrun.h ../src/json.h $(wildcard j*.inc) thm.inc pre.inc ../src/crimlab.h ../DESCRIPTION | $(B)
	$(CXX) $(WARN) $(CPPFLAGS) $(CXXFLAGS) $(JOPT) -fPIC -I../src -c $< -o $@

# A statically linked program cannot load a shared library (with glibc it
# stops at dlopen): --update then asks through the curl program alone.
ifneq ($(findstring -static,$(LDFLAGS)),)
UPDDEF = -DCLAB_NODL
endif

$(B)/main.o: main.cpp jrun.h ../src/json.h ../src/crimlab.h ../DESCRIPTION | $(B)
	$(CXX) $(WARN) $(CPPFLAGS) $(CXXFLAGS) $(UPDDEF) -I../src -DCLAB_VER='"$(VER)"' -c $< -o $@

# the proof table of the R binding
thm.inc: ../src/bind.cpp thm.sh
	sh thm.sh > $@

# the preset table of the R package
pre.inc: ../inst/presets/presets.tsv pre.sh
	sh pre.sh > $@

# the publishers' codes of the presets (prechk())
codes.inc: ../inst/presets/codes.tsv codes.sh
	sh codes.sh > $@
$(B)/jrun.o: codes.inc
# the version yoyhtml() writes in its page
$(B)/jrun.o: CPPFLAGS += -DCLAB_VER='"$(VER)"'

$(B):
	mkdir -p $@

libcrimlab.a: $(KO)
	rm -f $@
	$(AR) rcs $@ $(KO)

$(SO): $(KO)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOFLAGS) -o $@ $(KO) $(LDLIBS)

crimlab: $(B)/main.o libcrimlab.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(B)/main.o libcrimlab.a $(LDLIBS)

check: crimlab
	sh test.sh ./crimlab

install: all
	install -d $(DESTDIR)$(BINDIR)
	install -d $(DESTDIR)$(LIBDIR)
	install -d $(DESTDIR)$(INCDIR)/crimlab
	install -d $(DESTDIR)$(MANDIR)/man1
	install -m 755 crimlab $(DESTDIR)$(BINDIR)/crimlab
	install -m 755 $(SO) $(DESTDIR)$(LIBDIR)/$(SO)
	ln -sf $(SO) $(DESTDIR)$(LIBDIR)/$(if $(filter Darwin,$(shell uname -s)),libcrimlab.dylib,libcrimlab.so)
	install -m 644 ../src/crimlab.h jrun.h $(DESTDIR)$(INCDIR)/crimlab/
	install -m 644 crimlab.1 $(DESTDIR)$(MANDIR)/man1/crimlab.1

clean:
	rm -rf $(B) crimlab libcrimlab.a $(SO)

.PHONY: all check install clean
