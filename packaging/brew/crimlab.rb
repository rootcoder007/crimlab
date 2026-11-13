# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Homebrew formula (macOS and Linux) for the tap rootcoder007/crimlab
# (the repository rootcoder007/homebrew-crimlab):
#   brew install rootcoder007/crimlab/crimlab
# Builds the crimlab command from the release's source tarball. tap.yml
# sets url and sha256 for each published release; `brew install --HEAD`
# builds main.
class Crimlab < Formula
  desc "Criminology and sociolegal statistics from the command-line"
  homepage "https://github.com/rootcoder007/crimlab"
  url "https://github.com/rootcoder007/crimlab/releases/download/v0.0.0/crimlab-0.0.0.tar.gz"
  sha256 "0000000000000000000000000000000000000000000000000000000000000000"
  license "Apache-2.0" => { with: "LLVM-exception" }
  head "https://github.com/rootcoder007/crimlab.git", branch: "main"

  def install
    system "make", "-C", "cli", "crimlab", "CXX=#{ENV.cxx}", "CXXFLAGS=-O2", "WERROR="
    system "make", "-C", "cli", "check"
    bin.install "cli/crimlab"
    man1.install "cli/crimlab.1"
  end

  test do
    assert_match "0.1475205", shell_output("#{bin}/crimlab orrisk --or 2 --rate 0.1 --share 0.3")
    assert_match "\"nfail\":0", shell_output("#{bin}/crimlab selftest -j")
  end
end
