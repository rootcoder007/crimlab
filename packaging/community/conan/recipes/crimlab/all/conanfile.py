import os

from conan import ConanFile
from conan.errors import ConanInvalidConfiguration
from conan.tools.apple import is_apple_os
from conan.tools.build import check_min_cppstd
from conan.tools.files import chdir, copy, get
from conan.tools.gnu import Autotools, AutotoolsToolchain
from conan.tools.layout import basic_layout

required_conan_version = ">=2.1"


class CrimlabConan(ConanFile):
    name = "crimlab"
    description = "Statistical toolkit for criminology and sociolegal research"
    license = "Apache-2.0 WITH LLVM-exception"
    url = "https://github.com/conan-io/conan-center-index"
    homepage = "https://crimlab.dev/"
    topics = ("statistics", "criminology", "forensic-science", "causal-inference", "sociolegal-research")
    package_type = "library"
    settings = "os", "arch", "compiler", "build_type"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
    }

    @property
    def _soname(self):
        return "libcrimlab.0.dylib" if is_apple_os(self) else "libcrimlab.so.0"

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        basic_layout(self, src_folder="src")

    def validate(self):
        check_min_cppstd(self, 17)
        if self.settings.os == "Windows":
            raise ConanInvalidConfiguration(
                "crimlab builds with GNU make; its Makefile has no Windows build")

    def source(self):
        get(self, **self.conan_data["sources"][self.version], strip_root=True)

    def generate(self):
        tc = AutotoolsToolchain(self)
        tc.generate()

    def build(self):
        args = ["WERROR="]
        if is_apple_os(self):
            args.append(f"SOFLAGS=-dynamiclib -install_name @rpath/{self._soname}")
        target = self._soname if self.options.shared else "libcrimlab.a"
        autotools = Autotools(self)
        with chdir(self, os.path.join(self.source_folder, "cli")):
            autotools.make(target=target, args=args)

    def package(self):
        for f in ("LICENSE", "NOTICE"):
            copy(self, f, self.source_folder, os.path.join(self.package_folder, "licenses"))
        inc = os.path.join(self.package_folder, "include", "crimlab")
        copy(self, "crimlab.h", os.path.join(self.source_folder, "src"), inc)
        copy(self, "jrun.h", os.path.join(self.source_folder, "cli"), inc)
        lib = os.path.join(self.package_folder, "lib")
        cli = os.path.join(self.source_folder, "cli")
        if self.options.shared:
            copy(self, self._soname, cli, lib)
            link = "libcrimlab.dylib" if is_apple_os(self) else "libcrimlab.so"
            os.symlink(self._soname, os.path.join(lib, link))
        else:
            copy(self, "libcrimlab.a", cli, lib)

    def package_info(self):
        self.cpp_info.libs = ["crimlab"]
        if self.settings.os in ("Linux", "FreeBSD"):
            self.cpp_info.system_libs = ["m", "pthread", "dl"]
