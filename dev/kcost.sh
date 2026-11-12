#!/bin/sh
# An instrumented build for dev/cost.R: a copy of the package in which each
# kernel call of the binding (every kchk(clab_...) of src/bind.cpp) is timed
# with std::chrono::steady_clock, the seconds summed in a counter that
# .Call(crimlab:::C_clab_kt, reset) returns. The
# shipped sources are not changed.
#   sh dev/kcost.sh LIB        (from the package root; installs into LIB)
set -eu
lib=$1
t=$(mktemp -d)
trap 'rm -rf "$t"' EXIT
mkdir -p "$lib" "$t/crimlab"
tar cf - --exclude=.git --exclude='*.o' --exclude='*.so' . | tar xf - -C "$t/crimlab"
s=$t/crimlab/src
# the macro sits after the using-declaration, so rchk.h keeps the function
awk '{ print }
/^using clab::kchk;$/ {
  print "#include <chrono>"
  print "double clab_kt = 0;"
  print "#define kchk(c_, s_) do { const auto t0_ = std::chrono::steady_clock::now(); \\"
  print "  const int r_ = (c_); clab_kt += std::chrono::duration<double>( \\"
  print "  std::chrono::steady_clock::now() - t0_).count(); clab::kchk(r_, s_); } while (0)"
  print "extern \"C\" SEXP C_clab_kt(SEXP r) {"
  print "  SEXP v = Rf_ScalarReal(clab_kt); if (Rf_asLogical(r) == TRUE) clab_kt = 0; return v; }"
}' "$s/bind.cpp" > "$s/bind.tmp" && mv "$s/bind.tmp" "$s/bind.cpp"
awk '/^SEXP C_clab_tlen\(SEXP, SEXP\);$/ { print; print "SEXP C_clab_kt(SEXP);"; next }
/^    \{nullptr, nullptr, 0\}\};$/ { print "    {\"C_clab_kt\", fn(&C_clab_kt), 1},"; print; next }
{ print }' "$s/glue.cpp" > "$s/glue.tmp" && mv "$s/glue.tmp" "$s/glue.cpp"
grep -q C_clab_kt "$s/bind.cpp" && grep -q '"C_clab_kt"' "$s/glue.cpp"
R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null nice -n 10 \
  R CMD INSTALL --no-test-load --library="$lib" "$t/crimlab" >"$t/log" 2>&1 ||
  { tail -30 "$t/log"; exit 1; }
echo "kcost: instrumented crimlab in $lib"
