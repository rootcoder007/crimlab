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
WERROR ?= -Werror
WARN = -std=c++17 -Wall -Wextra -pedantic $(WERROR)
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include
MANDIR ?= $(PREFIX)/share/man
SOVER = 0
VER := $(shell sed -n 's/^Version: *//p' ../DESCRIPTION)
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

all: crimlab libcrimlab.a $(SO)

$(B)/%.o: ../src/%.cpp ../src/*.h | $(B)
	$(CXX) $(WARN) $(CPPFLAGS) $(CXXFLAGS) -fPIC -c $< -o $@

$(B)/jrun.o: jrun.cpp jrun.h json.h $(wildcard j*.inc) thm.inc ../src/crimlab.h | $(B)
	$(CXX) $(WARN) $(CPPFLAGS) $(CXXFLAGS) -fPIC -I../src -c $< -o $@

$(B)/main.o: main.cpp jrun.h json.h ../src/crimlab.h ../DESCRIPTION | $(B)
	$(CXX) $(WARN) $(CPPFLAGS) $(CXXFLAGS) -I../src -DCLAB_VER='"$(VER)"' -c $< -o $@

# the proof table of the R binding
thm.inc: ../src/bind.cpp thm.sh
	sh thm.sh > $@

$(B):
	mkdir -p $@

libcrimlab.a: $(KO)
	rm -f $@
	$(AR) rcs $@ $(KO)

$(SO): $(KO)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOFLAGS) -o $@ $(KO)

crimlab: $(B)/main.o libcrimlab.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(B)/main.o libcrimlab.a

check: crimlab
	sh test.sh ./crimlab

install: all
	install -d $(DESTDIR)$(BINDIR) $(DESTDIR)$(LIBDIR) $(DESTDIR)$(INCDIR)/crimlab $(DESTDIR)$(MANDIR)/man1
	install -m 755 crimlab $(DESTDIR)$(BINDIR)/crimlab
	install -m 755 $(SO) $(DESTDIR)$(LIBDIR)/$(SO)
	ln -sf $(SO) $(DESTDIR)$(LIBDIR)/$(if $(filter Darwin,$(shell uname -s)),libcrimlab.dylib,libcrimlab.so)
	install -m 644 ../src/crimlab.h jrun.h $(DESTDIR)$(INCDIR)/crimlab/
	install -m 644 crimlab.1 $(DESTDIR)$(MANDIR)/man1/crimlab.1

clean:
	rm -rf $(B) crimlab libcrimlab.a $(SO)

.PHONY: all check install clean
