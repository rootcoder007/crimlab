# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Fedora and EPEL RPM (COPR builds it from the repository through
# .copr/Makefile, which sets Version from DESCRIPTION).
Name:           crimlab
Version:        0.0.0.9000
Release:        1%{?dist}
Summary:        Criminology and sociolegal statistics: library and command
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
%make_build -C cli WERROR= PREFIX=%{_prefix} LIBDIR=%{_libdir}

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
* Thu Nov 12 2026 rootcoder007 <278967282+rootcoder007@users.noreply.github.com> - 0.0.0.9000-1
- First package
