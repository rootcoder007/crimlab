// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The byte checks of the exposure file of vrfpol() (R/vrfp.R): the end of
// a bzip2 stream, the end of a gzip member, and the encoding a text file
// is in. The file is read and decompressed by the caller.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::int64_t;

// bit i (from 0, most significant first) of the bytes t
int bit(const unsigned char *t, int64_t i) {
  return (t[i / 8] >> (7 - static_cast<int>(i % 8))) & 1;
}

std::uint32_t le32(const unsigned char *t) {
  return static_cast<std::uint32_t>(t[0]) | static_cast<std::uint32_t>(t[1]) << 8 |
         static_cast<std::uint32_t>(t[2]) << 16 | static_cast<std::uint32_t>(t[3]) << 24;
}

}  // namespace

extern "C" {

int clab_bz2end(const unsigned char *t, int64_t n) CLAB_TRY {
  if (n != 11) return 0;
  static const unsigned char mg[] = {0x17, 0x72, 0x45, 0x38, 0x50, 0x90};
  for (int64_t p = 1; p <= 8; ++p) {
    bool ok = true;
    for (int64_t k = 0; k < 48 && ok; ++k) ok = bit(t, p + k) == bit(mg, k);
    for (int64_t k = p + 80; k < 88 && ok; ++k) ok = bit(t, k) == 0;
    if (ok) return 1;
  }
  return 0;
} CLAB_CATCH(0)

int clab_gzend(const unsigned char *d, int64_t n, const unsigned char *t, int64_t nt) CLAB_TRY {
  if (nt != 8 || n < 0) return 0;
  const std::uint32_t len = le32(t + 4);
  if (static_cast<int64_t>(len) > n) return 0;
  const std::uint32_t crc = clab_crc32(d + (n - static_cast<int64_t>(len)), len, 0);
  return crc == le32(t) ? 1 : 0;
} CLAB_CATCH(0)

int clab_txtenc(const unsigned char *r, int64_t n, int *bom) CLAB_TRY {
  *bom = 0;
  static const unsigned char bm[4][4] = {
      {0xff, 0xfe, 0, 0}, {0, 0, 0xfe, 0xff}, {0xff, 0xfe, 0, 0}, {0xfe, 0xff, 0, 0}};
  static const int bl[4] = {4, 4, 2, 2};
  for (int e = 0; e < 4; ++e) {
    bool ok = n >= bl[e];
    for (int k = 0; k < bl[e] && ok; ++k) ok = r[k] == bm[e][k];
    if (ok) {
      *bom = bl[e];
      return e + 1;
    }
  }
  bool nul = false;
  for (int64_t i = 0; i < n && !nul; ++i) nul = r[i] == 0;
  if (!nul) return 0;
  // the NUL bytes of the first 2000 on the odd positions (from 1) only:
  // UTF-16BE; on the even ones: UTF-16LE
  const int64_t m = n < 2000 ? n : 2000;
  if (m < 2) return -1;
  double z0 = 0, z1 = 0, n0 = 0, n1 = 0;
  for (int64_t i = 0; i < m; ++i) {
    if (i % 2) {
      n0 += 1;
      z0 += r[i] == 0 ? 1 : 0;
    } else {
      n1 += 1;
      z1 += r[i] == 0 ? 1 : 0;
    }
  }
  z0 /= n0;
  z1 /= n1;
  if (z0 > 0.9 && z1 < 0.1) return 3;
  if (z1 > 0.9 && z0 < 0.1) return 4;
  return -1;
} CLAB_CATCH(-2)

}  // extern "C"
