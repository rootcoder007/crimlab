// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Text to plain ASCII (R/txt.R): bytes that are not UTF-8 dropped, the
// common Latin accents and ligatures folded by a fixed table, every other
// character outside the printable ASCII range (space to tilde) dropped.
// The same bytes on every system, whatever its locale.
#include <cstdint>
#include <cstring>

#include "crimlab.h"
#include "kutil.h"

namespace {

// One code point from s[i..n) into *cp; the bytes it takes, 0 when s[i]
// does not start a valid UTF-8 sequence (shortest form, no surrogates).
int dec(const unsigned char *s, int i, int n, std::uint32_t *cp) {
  const unsigned c = s[i];
  if (c < 0x80) {
    *cp = c;
    return 1;
  }
  int k = 0;
  std::uint32_t v = 0, lo = 0;
  if (c >= 0xC2 && c <= 0xDF) {
    k = 2;
    v = c & 0x1F;
    lo = 0x80;
  } else if (c >= 0xE0 && c <= 0xEF) {
    k = 3;
    v = c & 0x0F;
    lo = 0x800;
  } else if (c >= 0xF0 && c <= 0xF4) {
    k = 4;
    v = c & 0x07;
    lo = 0x10000;
  } else {
    return 0;
  }
  if (i + k > n) return 0;
  for (int j = 1; j < k; ++j) {
    const unsigned b = s[i + j];
    if ((b & 0xC0) != 0x80) return 0;
    v = (v << 6) | (b & 0x3F);
  }
  if (v < lo || v > 0x10FFFF || (v >= 0xD800 && v <= 0xDFFF)) return 0;
  *cp = v;
  return k;
}

// The folding of one code point into o: the number of ASCII bytes, 0 when
// it has none.
int fold(std::uint32_t c, char *o) {
  static const char *const l1 =  // U+00C0 to U+00FF
      "AAAAAA_CEEEEIIIIDNOOOOO_OUUUUY__aaaaaa_ceeeeiiiidnooooo_ouuuuy_y";
  static const std::uint32_t two[] = {0xDF, 0xC6, 0xE6, 0xDE, 0xFE, 0x152, 0x153};
  static const char *const tt[] = {"ss", "AE", "ae", "Th", "th", "OE", "oe"};
  for (int k = 0; k < 7; ++k) {
    if (two[k] == c) {
      o[0] = tt[k][0];
      o[1] = tt[k][1];
      return 2;
    }
  }
  if (c >= 0xC0 && c <= 0xFF) {
    o[0] = l1[c - 0xC0];
    return o[0] == '_' ? 0 : 1;
  }
  static const std::uint32_t ext[] = {
      0x100, 0x101, 0x104, 0x105, 0x106, 0x107, 0x10C, 0x10D, 0x10E, 0x10F, 0x110, 0x111, 0x112, 0x113,
      0x118, 0x119, 0x11A, 0x11B, 0x11E, 0x11F, 0x130, 0x131, 0x141, 0x142, 0x143, 0x144, 0x147, 0x148,
      0x14C, 0x14D, 0x150, 0x151, 0x158, 0x159, 0x15A, 0x15B, 0x15E, 0x15F, 0x160, 0x161, 0x164, 0x165,
      0x16A, 0x16B, 0x16E, 0x16F, 0x170, 0x171, 0x179, 0x17A, 0x17B, 0x17C, 0x17D, 0x17E};
  static const char *const to = "AaAaCcCcDdDdEeEeEeGgIiLlNnNnOoOoRrSsSsSsTtUuUuUuZzZzZz";
  for (std::size_t k = 0; k < sizeof ext / sizeof ext[0]; ++k) {
    if (ext[k] == c) {
      o[0] = to[k];
      return 1;
    }
  }
  return 0;
}

}  // namespace

extern "C" {

int clab_utf8ok(const char *s, int n) CLAB_TRY {
  const unsigned char *u = reinterpret_cast<const unsigned char *>(s);
  std::uint32_t cp = 0;
  for (int i = 0; i < n;) {
    const int k = dec(u, i, n, &cp);
    if (k == 0) return 0;
    i += k;
  }
  return 1;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_ascii(const char *s, int n, char *out) CLAB_TRY {
  const unsigned char *u = reinterpret_cast<const unsigned char *>(s);
  int m = 0;
  std::uint32_t cp = 0;
  for (int i = 0; i < n;) {
    const int k = dec(u, i, n, &cp);
    if (k == 0) {
      ++i;
      continue;
    }
    i += k;
    if (cp < 0x80) {
      if (cp >= 32 && cp <= 126) out[m++] = static_cast<char>(cp);
      continue;
    }
    m += fold(cp, out + m);
  }
  return m;
} CLAB_CATCH(CLAB_ENOMEM)

}  // extern "C"
