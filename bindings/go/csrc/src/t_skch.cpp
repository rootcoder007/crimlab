// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/skch.cpp: SHA-256 against the FIPS 180-4 examples, moments
// against their definition, the reservoir and the sketch.
#include <cmath>
#include <cstring>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_skch(tk &T, int) {
  // FIPS 180-4, appendix: "abc" and the 56-byte two-block message
  const unsigned char want1[4] = {0xba, 0x78, 0x16, 0xbf}, want2[4] = {0x24, 0x8d, 0x6a, 0x61};
  unsigned char h[32];
  clab_sha256(reinterpret_cast<const unsigned char *>("abc"), 3, h);
  tk_true(T, "skch.sha1", std::memcmp(h, want1, 4) == 0 && h[31] == 0xad);
  const char *m2 = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
  clab_sha256(reinterpret_cast<const unsigned char *>(m2), std::strlen(m2), h);
  tk_true(T, "skch.sha2", std::memcmp(h, want2, 4) == 0 && h[31] == 0xc1);
  // CRC-32 check value of "123456789" (RFC 1952), in one call and in two
  const unsigned char *c9 = reinterpret_cast<const unsigned char *>("123456789");
  tk_true(T, "skch.crc", clab_crc32(c9, 9, 0) == 0xCBF43926u && clab_crc32(c9, 0, 0) == 0);
  tk_true(T, "skch.crc2", clab_crc32(c9 + 4, 5, clab_crc32(c9, 4, 0)) == 0xCBF43926u);
  // moments of 1, 2, 3, 4 in two blocks equal those of all four
  const double x[4] = {1, 2, 3, 4};
  double a[5], b[5], m[5], all[5], s[6];
  clab_momacc(x, 2, a);
  clab_momacc(x + 2, 2, b);
  clab_mommrg(a, b, m);
  clab_momacc(x, 4, all);
  tk_veq(T, "skch.merge", m, all, 5, 1e-15);
  clab_momstat(all, s);
  tk_eq(T, "skch.var", s[2], 5.0 / 3, 1e-15);
  tk_eq(T, "skch.skew", s[4], 0, 1e-15);
  clab_momacc(x, 0, a);
  tk_eq(T, "skch.empty", a[0], 0, 0);
  // the reservoir of a stream no longer than k is the stream
  double r[3];
  clab_rsv(3, 5, 1, r);
  tk_eq(T, "skch.rsv", r[2], 3, 0);
  tk_ieq(T, "skch.rsvk", clab_rsv(3, 2e8, 1, r), CLAB_EARG);
  // one value fills one register; the estimate of one value is near 1
  std::vector<int> reg(16, 0);
  clab_hlladd("a", 1, 4, reg.data());
  tk_near(T, "skch.hll", clab_hllcnt(reg.data(), 16), 1, 0.05);
  tk_ieq(T, "skch.hllp", clab_hlladd("a", 1, 3, reg.data()), CLAB_EARG);
  tk_eq(T, "skch.hll0", clab_hllcnt(reg.data(), 0), 0, 0);
}
