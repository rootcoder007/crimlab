// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One-pass summaries of data larger than memory: moment accumulators merged
// exactly (Chan, Golub and LeVeque 1983; Terriberry's third and fourth
// moments), reservoir sampling (Vitter 1985, Algorithm R) with Lemire's
// (2019) unbiased bounded draw on std::mt19937_64, and the HyperLogLog
// distinct-value sketch (Flajolet, Fusy, Gandouet and Meunier 2007) on
// SHA-256 (FIPS 180-4), and the CRC-32 of gzip (RFC 1952).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::uint32_t;
using std::uint64_t;

const uint32_t K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

void block(uint32_t *h, const unsigned char *p) {
  uint32_t w[64];
  for (int i = 0; i < 16; ++i) {
    w[i] = static_cast<uint32_t>(p[4 * i]) << 24 | static_cast<uint32_t>(p[4 * i + 1]) << 16 |
           static_cast<uint32_t>(p[4 * i + 2]) << 8 | static_cast<uint32_t>(p[4 * i + 3]);
  }
  for (int i = 16; i < 64; ++i) {
    const uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
    const uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
    w[i] = w[i - 16] + s0 + w[i - 7] + s1;
  }
  uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], k = h[7];
  for (int i = 0; i < 64; ++i) {
    const uint32_t t1 = k + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + K[i] + w[i];
    const uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
    k = g;
    g = f;
    f = e;
    e = d + t1;
    d = c;
    c = b;
    b = a;
    a = t1 + t2;
  }
  h[0] += a;
  h[1] += b;
  h[2] += c;
  h[3] += d;
  h[4] += e;
  h[5] += f;
  h[6] += g;
  h[7] += k;
}

// Lemire's (2019) unbiased draw from [0, s) on 64-bit outputs.
uint64_t bdraw(std::mt19937_64 &g, uint64_t s) {
  auto mul = [](uint64_t x, uint64_t y, uint64_t &hi, uint64_t &lo) {
    const uint64_t x0 = x & 0xFFFFFFFFULL, x1 = x >> 32, y0 = y & 0xFFFFFFFFULL, y1 = y >> 32;
    const uint64_t p00 = x0 * y0, p01 = x0 * y1, p10 = x1 * y0, p11 = x1 * y1;
    const uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFFULL) + (p10 & 0xFFFFFFFFULL);
    hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
    lo = (mid << 32) | (p00 & 0xFFFFFFFFULL);
  };
  uint64_t hi, lo;
  mul(g(), s, hi, lo);
  if (lo < s) {
    const uint64_t t = (0 - s) % s;
    while (lo < t) mul(g(), s, hi, lo);
  }
  return hi;
}

}  // namespace

extern "C" {

void clab_sha256(const unsigned char *d, std::size_t n, unsigned char *out) {
  uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  std::size_t i = 0;
  for (; i + 64 <= n; i += 64) block(h, d + i);
  unsigned char t[128] = {0};
  const std::size_t r = n - i;
  if (r > 0) std::memcpy(t, d + i, r);
  t[r] = 0x80;
  const std::size_t tl = r + 9 <= 64 ? 64 : 128;
  const uint64_t bits = static_cast<uint64_t>(n) * 8;
  for (int k = 0; k < 8; ++k) t[tl - 1 - k] = static_cast<unsigned char>(bits >> (8 * k));
  block(h, t);
  if (tl == 128) block(h, t + 64);
  for (int k = 0; k < 8; ++k) {
    for (int b = 0; b < 4; ++b) out[4 * k + b] = static_cast<unsigned char>(h[k] >> (24 - 8 * b));
  }
}

uint32_t clab_crc32(const unsigned char *d, std::size_t n, uint32_t crc) {
  // the reflected polynomial 0xEDB88320, one byte at a time
  crc = ~crc;
  for (std::size_t i = 0; i < n; ++i) {
    crc ^= d[i];
    for (int k = 0; k < 8; ++k) crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
  }
  return ~crc;
}

int clab_momacc(const double *x, double n, double *out) {
  for (int i = 0; i < 5; ++i) out[i] = 0;
  if (n < 1) return CLAB_OK;
  const std::size_t m = static_cast<std::size_t>(n);
  double mean = 0;
  for (std::size_t i = 0; i < m; ++i) mean += (x[i] - mean) / (static_cast<double>(i) + 1.0);
  double r = 0;
  for (std::size_t i = 0; i < m; ++i) r += x[i] - mean;
  mean += r / n;
  double m2 = 0, m3 = 0, m4 = 0;
  for (std::size_t i = 0; i < m; ++i) {
    const double d = x[i] - mean, d2 = d * d;
    m2 += d2;
    m3 += d2 * d;
    m4 += d2 * d2;
  }
  out[0] = n;
  out[1] = mean;
  out[2] = m2;
  out[3] = m3;
  out[4] = m4;
  return CLAB_OK;
}

int clab_mommrg(const double *a, const double *b, double *out) {
  const double na = a[0], nb = b[0];
  if (na <= 0 || nb <= 0) {
    const double *s = na <= 0 ? b : a;
    for (int i = 0; i < 5; ++i) out[i] = s[i];
    return CLAB_OK;
  }
  const double n = na + nb, dl = b[1] - a[1], d2 = dl * dl, d3 = d2 * dl, d4 = d2 * d2;
  const double mean = a[1] + nb * dl / n;
  const double M2 = a[2] + b[2] + d2 * na * nb / n;
  const double M3 = a[3] + b[3] + d3 * na * nb * (na - nb) / (n * n) + 3.0 * dl * (na * b[2] - nb * a[2]) / n;
  const double M4 = a[4] + b[4] + d4 * na * nb * (na * na - na * nb + nb * nb) / (n * n * n) +
                    6.0 * d2 * (na * na * b[2] + nb * nb * a[2]) / (n * n) + 4.0 * dl * (na * b[3] - nb * a[3]) / n;
  out[0] = n;
  out[1] = mean;
  out[2] = M2;
  out[3] = M3;
  out[4] = M4;
  return CLAB_OK;
}

int clab_rsv(double n, double k, double seed, double *out) {
  if (!(n >= 0) || !(k >= 0) || !(seed >= 0) || n > 9007199254740992.0 || k > 1e8 || seed >= 9007199254740992.0) {
    return CLAB_EARG;
  }
  const uint64_t nn = static_cast<uint64_t>(n);
  uint64_t kk = static_cast<uint64_t>(k);
  if (kk > nn) kk = nn;
  std::mt19937_64 g(static_cast<unsigned long long>(seed));
  std::vector<uint64_t> res(static_cast<std::size_t>(kk));
  for (uint64_t i = 0; i < kk; ++i) res[static_cast<std::size_t>(i)] = i;
  for (uint64_t j = kk; j < nn; ++j) {
    if ((j & 0xFFFFF) == 0 && kutil::stop(1048576.0)) return CLAB_EINTR;
    const uint64_t t = bdraw(g, j + 1);
    if (t < kk) res[static_cast<std::size_t>(t)] = j;
  }
  std::sort(res.begin(), res.end());
  for (uint64_t i = 0; i < kk; ++i) out[i] = static_cast<double>(res[static_cast<std::size_t>(i)] + 1);
  return CLAB_OK;
}

int clab_hlladd(const char *s, std::size_t len, int p, int *reg) {
  if (p < 4 || p > 20) return CLAB_EARG;
  unsigned char h[32];
  clab_sha256(reinterpret_cast<const unsigned char *>(s), len, h);
  uint64_t v = 0;
  for (int i = 0; i < 8; ++i) v = (v << 8) | h[i];
  const std::size_t idx = static_cast<std::size_t>(v >> (64 - p));
  const uint64_t rest = (v << p) | (static_cast<uint64_t>(1) << (p - 1));
  int rank = 1;
  uint64_t mask = static_cast<uint64_t>(1) << 63;
  while (rank <= 64 - p && (rest & mask) == 0) {
    ++rank;
    mask >>= 1;
  }
  if (rank > reg[idx]) reg[idx] = rank;
  return CLAB_OK;
}

double clab_hllcnt(const int *r, std::size_t m) {
  if (m == 0) return 0;
  double sum = 0;
  std::size_t zeros = 0;
  for (std::size_t i = 0; i < m; ++i) {
    sum += std::pow(2.0, -static_cast<double>(r[i]));
    if (r[i] == 0) ++zeros;
  }
  const double md = static_cast<double>(m);
  const double al = m == 16 ? 0.673 : m == 32 ? 0.697 : m == 64 ? 0.709 : 0.7213 / (1.0 + 1.079 / md);
  double est = al * md * md / sum;
  if (est <= 2.5 * md && zeros > 0) est = md * std::log(md / static_cast<double>(zeros));
  return est;
}

}  // extern "C"

extern "C" void clab_momstat(const double *s, double *o) {
  const double n = s[0], nan = kutil::nan();
  const double m2 = n > 0 ? s[2] / n : nan;
  o[0] = n;
  o[1] = n > 0 ? s[1] : nan;
  o[2] = n > 1 ? s[2] / (n - 1) : nan;
  o[3] = std::sqrt(o[2]);
  o[4] = n > 2 && s[2] > 0 ? (s[3] / n) / std::pow(m2, 1.5) : nan;
  o[5] = n > 3 && s[2] > 0 ? (s[4] / n) / (m2 * m2) - 3 : nan;
}
