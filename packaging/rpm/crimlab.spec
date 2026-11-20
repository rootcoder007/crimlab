# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Fedora and EPEL RPM (COPR builds it from the repository through
# .copr/Makefile, which sets Version from DESCRIPTION).
Name:           crimlab
Version:        0.0.4
Release:        1%{?dist}
Summary:        Statistical Toolkit for Criminology and Sociolegal Research: library and command
License:        Apache-2.0 WITH LLVM-exception
URL:            https://github.com/rootcoder007/crimlab
Source0:        %{url}/archive/refs/tags/v%{version}/%{name}-%{version}.tar.gz

BuildRequires:  gcc-c++
BuildRequires:  make

%description
Statistical methods for courts, corrections, policing, forensic science
and crime places, the compiled kernels of the R and Python packages
crimlab. This package has the shared library libcrimlab and the crimlab
command, which runs every function of the package on JSON or CSV input.

%package        devel
Summary:        Headers for %{name}
Requires:       %{name}%{?_isa} = %{version}-%{release}

%description    devel
The C headers of libcrimlab: the kernels (crimlab.h) and the
call-by-name interface with JSON in and out (jrun.h).

%prep
%autosetup

%build
%set_build_flags
# The riscv64 builders are emulated, and the JSON layer is one large file:
# at -O2 with full debug information it alone took them over an hour. It
# is built at -O1 with line tables there; the kernels keep the flags of
# the distribution.
%ifarch riscv64
%global jopt JOPT="-O1 -g1"
%endif
%make_build -C cli WERROR= PREFIX=%{_prefix} LIBDIR=%{_libdir} %{?jopt}

%install
%make_install -C cli WERROR= PREFIX=%{_prefix} LIBDIR=%{_libdir} MANDIR=%{_mandir}

%check
make -C cli check WERROR=

%files
%license LICENSE NOTICE
%doc README.md
%{_bindir}/crimlab
%{_libdir}/libcrimlab.so.0
%{_mandir}/man1/crimlab.1*

%files devel
%{_includedir}/crimlab/
%{_libdir}/libcrimlab.so

%changelog
* Fri Nov 20 2026 rootcoder007 <278967282+rootcoder007@users.noreply.github.com> - 0.0.4-1
- The loader and thread libraries are named at the link, for C libraries that keep them outside libc
- make install creates each directory in a command of its own
- No change to the functions

* Thu Nov 19 2026 rootcoder007 <278967282+rootcoder007@users.noreply.github.com> - 0.0.3-1
- The JSON layer is compiled at -O1 with line tables on riscv64, where it took the builders over an hour
- No change to the functions

* Wed Nov 18 2026 rootcoder007 <278967282+rootcoder007@users.noreply.github.com> - 0.0.2-1
- Every package is described by the title of the R package and says how to cite it
- New: plotsvg (the plots as SVG), plotcoh, rvwtim, entcmp, dmvar; the argument that names a method is mthd
- Fixed: a missing value stays missing on riscv64; 32-bit x86 is built with the arithmetic of SSE2

* Thu Nov 12 2026 rootcoder007 <278967282+rootcoder007@users.noreply.github.com> - 0.0.1-1
- First package
