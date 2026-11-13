// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Shortest round-trip decimal text for a double, without the C library:
// the fewest significant digits that read back as the same double (when
// two such strings are equally short, the one nearer the double), by the
// free-format algorithm of Burger and Dybvig (1996,
// doi:10.1145/231379.231397) in exact integer arithmetic. The same bytes
// on every system, whatever its printf, locale or rounding of products.
// Header only; nothing here allocates.
//
// Format: as Python's repr() without its ".0": fixed notation when the
// decimal exponent is from -4 to 15 (0.0001, 1234.5, 1e+16 is not),
// otherwise d.ddde+XX with at least two exponent digits; "NaN", "Inf",
// "-Inf"; "0" and "-0".
#ifndef CLAB_DFMT_H
#define CLAB_DFMT_H

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace dfmt {

// A non-negative integer of up to 48 words of 32 bits (1536 bits: the
// largest number below is about 2^1130).
struct Big {
  std::uint32_t w[48] = {};
  int n = 0;  // words in use; w[n - 1] != 0 unless n == 0
  void set(std::uint64_t v) noexcept {
    n = 0;
    while (v) {
      w[n++] = static_cast<std::uint32_t>(v);
      v >>= 32;
    }
  }
  void trim() noexcept {
    while (n > 0 && w[n - 1] == 0) --n;
  }
  void mul(std::uint32_t m) noexcept {
    std::uint64_t c = 0;
    for (int i = 0; i < n; ++i) {
      c += static_cast<std::uint64_t>(w[i]) * m;
      w[i] = static_cast<std::uint32_t>(c);
      c >>= 32;
    }
    if (c) w[n++] = static_cast<std::uint32_t>(c);
  }
  void shl(int b) noexcept {  // times 2^b
    if (n == 0) return;
    const int q = b / 32, r = b % 32;
    if (r) {
      std::uint32_t c = 0;
      for (int i = 0; i < n; ++i) {
        const std::uint32_t t = w[i] >> (32 - r);
        w[i] = (w[i] << r) | c;
        c = t;
      }
      if (c) w[n++] = c;
    }
    if (q) {
      for (int i = n - 1; i >= 0; --i) w[i + q] = w[i];
      for (int i = 0; i < q; ++i) w[i] = 0;
      n += q;
    }
  }
  void pow10(int k) noexcept {  // times 10^k
    for (; k >= 9; k -= 9) mul(1000000000U);
    for (; k > 0; --k) mul(10);
  }
  void sub(const Big &b) noexcept {  // b <= *this
    std::int64_t c = 0;
    for (int i = 0; i < n; ++i) {
      c += static_cast<std::int64_t>(w[i]) - (i < b.n ? b.w[i] : 0);
      w[i] = static_cast<std::uint32_t>(c);
      c = c < 0 ? -1 : 0;
    }
    trim();
  }
};

// Sign of (a + b) - c.
inline int cmpsum(const Big &a, const Big &b, const Big &c) noexcept {
  std::uint32_t s[49];
  const int m = a.n > b.n ? a.n : b.n;
  std::uint64_t k = 0;
  for (int i = 0; i < m; ++i) {
    k += static_cast<std::uint64_t>(i < a.n ? a.w[i] : 0) + (i < b.n ? b.w[i] : 0);
    s[i] = static_cast<std::uint32_t>(k);
    k >>= 32;
  }
  int n = m;
  if (k) s[n++] = static_cast<std::uint32_t>(k);
  if (n != c.n) return n < c.n ? -1 : 1;
  for (int i = n - 1; i >= 0; --i) {
    if (s[i] != c.w[i]) return s[i] < c.w[i] ? -1 : 1;
  }
  return 0;
}

inline int cmp(const Big &a, const Big &b) noexcept {
  if (a.n != b.n) return a.n < b.n ? -1 : 1;
  for (int i = a.n - 1; i >= 0; --i) {
    if (a.w[i] != b.w[i]) return a.w[i] < b.w[i] ? -1 : 1;
  }
  return 0;
}

// The text of x, NUL-terminated, in out (32 bytes is enough); returns
// its length.
inline int put(double x, char *out) noexcept {
  std::uint64_t u;
  std::memcpy(&u, &x, sizeof u);
  char *p = out;
  const bool neg = (u >> 63) != 0;
  const int be = static_cast<int>((u >> 52) & 0x7ff);
  std::uint64_t f = u & ((1ULL << 52) - 1);
  if (be == 0x7ff) {
    const char *t = f ? "NaN" : neg ? "-Inf" : "Inf";
    const std::size_t n = std::strlen(t);
    std::memcpy(out, t, n + 1);
    return static_cast<int>(n);
  }
  if (neg) *p++ = '-';
  if (be == 0 && f == 0) {
    *p++ = '0';
    *p = '\0';
    return static_cast<int>(p - out);
  }
  int e;
  if (be == 0) {
    e = -1074;
  } else {
    f |= 1ULL << 52;
    e = be - 1075;
  }
  // x = f 2^e. r / s is x scaled; mp and mm are half the gaps to the
  // neighbouring doubles above and below, on the same scale.
  const bool even = (f & 1) == 0;
  const bool edge = f == (1ULL << 52) && be > 1;  // the gap below is half
  Big r, s, mp, mm;
  r.set(f);
  s.set(1);
  mp.set(1);
  mm.set(1);
  if (e >= 0) {
    r.shl(e + (edge ? 2 : 1));
    s.shl(edge ? 2 : 1);
    mp.shl(e + (edge ? 1 : 0));
    mm.shl(e);
  } else {
    r.shl(edge ? 2 : 1);
    s.shl(edge ? 2 - e : 1 - e);
    if (edge) mp.shl(1);
  }
  // k, the decimal exponent: from the binary one, x >= 2^m, an estimate
  // of floor(m log10 2) (78913 / 2^18 is log10 2 to 3e-8) less one, so
  // never too large; the loop below raises it to the first k with the
  // high end of the interval below 10^k
  int m = e - 1;
  for (std::uint64_t t = f; t; t >>= 1) ++m;
  const long long q10 = m * 78913LL;
  int k = static_cast<int>(q10 >= 0 ? q10 >> 18 : -((-q10 + (1LL << 18) - 1) >> 18)) - 1;
  if (k >= 0) {
    s.pow10(k);
  } else {
    r.pow10(-k);
    mp.pow10(-k);
    mm.pow10(-k);
  }
  for (;;) {  // the high end of the interval must be below 10^k
    const int c = cmpsum(r, mp, s);
    if (even ? c < 0 : c <= 0) break;
    s.mul(10);
    ++k;
  }
  // digits
  char d[24];
  int nd = 0;
  for (;;) {
    r.mul(10);
    mp.mul(10);
    mm.mul(10);
    int q = 0;
    while (cmp(r, s) >= 0) {
      r.sub(s);
      ++q;
    }
    const int lo = cmp(r, mm), hi = cmpsum(r, mp, s);
    const bool t1 = even ? lo <= 0 : lo < 0;
    const bool t2 = even ? hi >= 0 : hi > 0;
    if (!t1 && !t2) {
      d[nd++] = static_cast<char>('0' + q);
      continue;
    }
    if (t1 && t2) {
      // both round-trip: the nearer to x, 2r against s
      const int h = cmpsum(r, r, s);
      if (h > 0 || (h == 0 && (q & 1))) ++q;
    } else if (t2) {
      ++q;
    }
    d[nd++] = static_cast<char>('0' + q);
    break;
  }
  // a final 10 carries
  for (int i = nd - 1; i > 0 && d[i] > '9'; --i) {
    d[i] = '0';
    ++d[i - 1];
  }
  if (d[0] > '9') {
    d[0] = '1';
    ++k;
  }
  while (nd > 1 && d[nd - 1] == '0') --nd;
  // x = 0.d * 10^k; its first digit has decimal exponent k - 1
  const int x10 = k - 1;
  if (x10 >= -4 && x10 < 16) {
    if (k <= 0) {
      *p++ = '0';
      *p++ = '.';
      for (int i = 0; i < -k; ++i) *p++ = '0';
      for (int i = 0; i < nd; ++i) *p++ = d[i];
    } else {
      for (int i = 0; i < k; ++i) *p++ = i < nd ? d[i] : '0';
      if (nd > k) {
        *p++ = '.';
        for (int i = k; i < nd; ++i) *p++ = d[i];
      }
    }
  } else {
    *p++ = d[0];
    if (nd > 1) {
      *p++ = '.';
      for (int i = 1; i < nd; ++i) *p++ = d[i];
    }
    *p++ = 'e';
    *p++ = x10 < 0 ? '-' : '+';
    const int a = x10 < 0 ? -x10 : x10;
    if (a >= 100) *p++ = static_cast<char>('0' + a / 100);
    *p++ = static_cast<char>('0' + a / 10 % 10);
    *p++ = static_cast<char>('0' + a % 10);
  }
  *p = '\0';
  return static_cast<int>(p - out);
}

}  // namespace dfmt

#endif
