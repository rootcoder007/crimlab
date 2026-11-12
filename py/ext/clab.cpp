// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The Python binding of libcrimlab, the module crimlab._clab. One entry
// point per exported function: each checks its arguments in the order and
// words of the R binding (src/bind.cpp, src/rchk.cpp), calls the kernel and
// builds the result from plain Python objects: a data frame is a dict of
// lists, a named vector a dict of floats, a list a dict. Names and labels
// are resolved in Python first (crimlab/_util.py): an entry point gets the
// values of a named argument already in the order of the groups, and the
// group names as a list (None when there are none).
//
// Numbers are read from NumPy's array interface without a copy when they
// are contiguous float64, and through the sequence protocol otherwise. (The
// stable ABI has the buffer protocol only from Python 3.11; this module is
// built for the stable ABI of 3.10.) Inside an entry point an error is a
// C++ exception; every Python object is held by an owner (Obj) that
// releases it on every path, and the entry point turns the exception into
// the Python one.
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "crimlab.h"

extern "C" int clab_selftest(char *sink, int cap, int *nfail);

namespace {

// The function an error belongs to, put in front of every message.
const char *who = "crimlab";

// Thrown once a Python exception is set.
struct Raised {};

// An owned reference.
class Obj {
 public:
  Obj() = default;
  explicit Obj(PyObject *p) : p_(p) {}
  ~Obj() { Py_XDECREF(p_); }
  Obj(Obj &&o) noexcept : p_(o.p_) { o.p_ = nullptr; }
  Obj &operator=(Obj &&o) noexcept {
    std::swap(p_, o.p_);
    return *this;
  }
  Obj(const Obj &) = delete;
  Obj &operator=(const Obj &) = delete;
  PyObject *get() const { return p_; }
  PyObject *release() {
    PyObject *t = p_;
    p_ = nullptr;
    return t;
  }

 private:
  PyObject *p_ = nullptr;
};

// A new reference, or the error it failed with.
Obj own(PyObject *p) {
  if (p == nullptr) throw Raised{};
  return Obj(p);
}

std::string vfmt(const char *f, std::va_list ap) {
  std::va_list aq;
  va_copy(aq, ap);
  std::vector<char> b(256);
  const int n = std::vsnprintf(b.data(), b.size(), f, ap);
  if (n >= static_cast<int>(b.size())) {
    b.resize(static_cast<std::size_t>(n) + 1);
    std::vsnprintf(b.data(), b.size(), f, aq);
  }
  va_end(aq);
  return std::string(b.data(), n > 0 ? static_cast<std::size_t>(n) : 0);
}

[[noreturn]] void raise(PyObject *type, const std::string &m) {
  PyErr_Format(type, "%s: %s", who, m.c_str());
  throw Raised{};
}

// ValueError (TypeError when type is true), "<who>: <message>".
[[noreturn]] void fail(bool type, const char *f, ...) {
  std::va_list ap;
  va_start(ap, f);
  const std::string m = vfmt(f, ap);
  va_end(ap);
  raise(type ? PyExc_TypeError : PyExc_ValueError, m);
}

void kchk(int code) {
  if (code != CLAB_OK) fail(false, "%s", clab_strerror(code));
}

// R's NA_real_: a NaN whose low word is 1954. Only the binding writes it,
// where R's result holds NA; it becomes None.
double mkna() {
  const std::uint64_t u = 0x7ff00000000007a2ULL;
  double d;
  std::memcpy(&d, &u, sizeof d);
  return d;
}
const double NA = mkna();
bool isna(double x) {
  std::uint64_t u;
  std::memcpy(&u, &x, sizeof u);
  return std::isnan(x) && (u & 0xffffffffULL) == 1954;
}
const double NaN = std::numeric_limits<double>::quiet_NaN();

// as.character() of one double, as R 4.3 and later writes it: the fewest
// significant digits (at most 15, or every digit of a whole number below
// 2^53) that give the same number, fixed notation unless scientific is
// shorter.
std::string rfmt(double x) {
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  if (x == 0) return "0";
  char b[64];
  const bool whole = x == std::floor(x) && std::fabs(x) <= 9007199254740992.0;
  const int top = whole ? 17 : 15;
  std::snprintf(b, sizeof b, "%.*e", top - 1, x);
  const double v = std::strtod(b, nullptr);
  int d = 1;
  for (; d < top; ++d) {
    std::snprintf(b, sizeof b, "%.*e", d - 1, x);
    if (std::strtod(b, nullptr) == v) break;
  }
  std::snprintf(b, sizeof b, "%.*e", d - 1, x);
  const int e = std::atoi(std::strchr(b, 'e') + 1);
  const int neg = x < 0;
  const int sci = neg + d + (d > 1) + (std::abs(e) >= 100 ? 5 : 4);
  const int rgt = std::max(0, d - e - 1);
  const int fix = neg + (e >= 0 ? e + 1 : 1) + (rgt > 0 ? rgt + 1 : 0);
  if (fix <= sci) std::snprintf(b, sizeof b, "%.*f", rgt, x);
  return b;
}

// The arithmetic helpers of the R binding (src/rchk.cpp), the same
// operations in the same order.
double rd(double x) {
  volatile double v = x;
  return v;
}

double nsum(const double *x, Py_ssize_t n) {
  double s = 0.0, c = 0.0;
  for (Py_ssize_t i = 0; i < n; ++i) {
    const double t = s + x[i];
    c += std::fabs(s) >= std::fabs(x[i]) ? (s - t) + x[i] : (x[i] - t) + s;
    s = t;
  }
  return std::isfinite(s) ? s + c : s;
}

double nmean(const double *x, Py_ssize_t n) {
  const double dn = static_cast<double>(n);
  double s = 0.0;
  for (Py_ssize_t i = 0; i < n; ++i) s += x[i];
  s /= dn;
  if (!std::isfinite(s)) return s;
  double t = 0.0;
  for (Py_ssize_t i = 0; i < n; ++i) t += x[i] - s;
  return s + t / dn;
}

// R's x ^ y for x = 10: R_pow() is pow() there.
double rpow(double x, double y) { return y == 2.0 ? x * x : std::pow(x, y); }

int ilen(double n, const char *nm) {
  if (n > INT_MAX) {
    if (nm == nullptr) fail(false, "%s", clab_strerror(CLAB_EARG));
    fail(false, "`%s` is too long for crimlab: more than %d values.", nm, INT_MAX);
  }
  return static_cast<int>(n);
}

// Reading arguments.

// What a reader found: numbers (or logicals), something else, or a matrix
// with more than one dimension above 1.
enum Got { NUM, NOTNUM, MATRIX };

// The array interface of an object: element kind and size, the number of
// elements, and the address and step of the one dimension above 1 (big
// counts such dimensions).
struct Ai {
  char kind = 0;
  int size = 0;
  int big = 0;
  Py_ssize_t n = 1, step = 0;
  const char *data = nullptr;
};

bool little() {
  const std::uint16_t x = 1;
  unsigned char c;
  std::memcpy(&c, &x, 1);
  return c == 1;
}

// False when o has no array interface this reader can use (none, no
// address, a byte order not the machine's): the caller then reads o as a
// sequence.
bool arrif(PyObject *o, Ai &a) {
  PyObject *raw = PyObject_GetAttrString(o, "__array_interface__");
  if (raw == nullptr) {
    PyErr_Clear();
    return false;
  }
  Obj ai(raw);
  if (!PyDict_Check(raw)) return false;
  PyObject *ts = PyDict_GetItemString(raw, "typestr");
  PyObject *sh = PyDict_GetItemString(raw, "shape");
  PyObject *da = PyDict_GetItemString(raw, "data");
  PyObject *st = PyDict_GetItemString(raw, "strides");
  if (ts == nullptr || !PyUnicode_Check(ts) || sh == nullptr || !PyTuple_Check(sh) ||
      da == nullptr || !PyTuple_Check(da) || PyTuple_Size(da) < 1) {
    return false;
  }
  const Py_ssize_t nd = PyTuple_Size(sh);
  if (st != nullptr && st != Py_None && (!PyTuple_Check(st) || PyTuple_Size(st) != nd)) {
    return false;
  }
  Py_ssize_t len = 0;
  const char *t = PyUnicode_AsUTF8AndSize(ts, &len);
  if (t == nullptr) throw Raised{};
  if (len < 2) return false;
  const bool native = t[0] == '|' || t[0] == '=' || (t[0] == '<') == little();
  if (!native) return false;
  a.kind = t[1];
  a.size = std::atoi(t + 2);
  a.step = a.size;
  for (Py_ssize_t i = 0; i < nd; ++i) {
    const Py_ssize_t d = PyLong_AsSsize_t(PyTuple_GetItem(sh, i));
    if (d == -1 && PyErr_Occurred()) throw Raised{};
    a.n *= d;
    if (d > 1) {
      ++a.big;
      if (st != nullptr && st != Py_None) {
        a.step = PyLong_AsSsize_t(PyTuple_GetItem(st, i));
        if (a.step == -1 && PyErr_Occurred()) throw Raised{};
      }
    }
  }
  a.data = static_cast<const char *>(PyLong_AsVoidPtr(PyTuple_GetItem(da, 0)));
  if (a.data == nullptr && PyErr_Occurred()) throw Raised{};
  return true;
}

// Element i of a numeric array as a double.
double aelt(const Ai &a, Py_ssize_t i) {
  const char *e = a.data + i * a.step;
#define CLAB_AS(T)            \
  {                           \
    T v;                      \
    std::memcpy(&v, e, sizeof v); \
    return static_cast<double>(v); \
  }
  if (a.kind == 'f') {
    if (a.size == 4) CLAB_AS(float)
    CLAB_AS(double)
  }
  if (a.kind == 'i') {
    if (a.size == 1) CLAB_AS(std::int8_t)
    if (a.size == 2) CLAB_AS(std::int16_t)
    if (a.size == 4) CLAB_AS(std::int32_t)
    CLAB_AS(std::int64_t)
  }
  if (a.size == 1) CLAB_AS(std::uint8_t)
  if (a.size == 2) CLAB_AS(std::uint16_t)
  if (a.size == 4) CLAB_AS(std::uint32_t)
  CLAB_AS(std::uint64_t)
#undef CLAB_AS
}

bool numkind(const Ai &a) {
  if (a.kind == 'f') return a.size == 4 || a.size == 8;
  return (a.kind == 'i' || a.kind == 'u') &&
         (a.size == 1 || a.size == 2 || a.size == 4 || a.size == 8);
}

// The name of an object's type.
std::string tname(PyObject *o) {
  Obj n = own(PyObject_GetAttrString(reinterpret_cast<PyObject *>(Py_TYPE(o)), "__name__"));
  const char *s = PyUnicode_AsUTF8AndSize(n.get(), nullptr);
  if (s == nullptr) throw Raised{};
  return s;
}

bool isbool(PyObject *o) {
  if (PyBool_Check(o)) return true;
  if (PyFloat_Check(o) || PyLong_Check(o)) return false;
  const std::string t = tname(o);
  return t == "bool_" || t == "bool";
}

// One element read as a number: None and pandas.NA are missing (NaN); a
// bool, text or anything float() refuses is not a number.
bool elt(PyObject *e, double &x) {
  if (e == Py_None) {
    x = NaN;
    return true;
  }
  if (PyUnicode_Check(e) || PyBytes_Check(e) || isbool(e)) return false;
  if (PyFloat_Check(e)) {
    x = PyFloat_AsDouble(e);
    return true;
  }
  if (PyLong_Check(e)) {
    x = PyLong_AsDouble(e);
    if (x == -1.0 && PyErr_Occurred()) {  // too large for a double
      PyErr_Clear();
      Obj z = own(PyLong_FromLong(0));
      const int neg = PyObject_RichCompareBool(e, z.get(), Py_LT);
      if (neg < 0) throw Raised{};
      x = neg ? -HUGE_VAL : HUGE_VAL;
    }
    return true;
  }
  if (tname(e) == "NAType") {
    x = NaN;
    return true;
  }
  PyObject *f = PyNumber_Float(e);
  if (f == nullptr) {
    PyErr_Clear();
    return false;
  }
  x = PyFloat_AsDouble(f);
  Py_DECREF(f);
  return true;
}

// Doubles read from a Python object: p points into the caller's array
// (no copy) or into v.
struct Dv {
  const double *p = nullptr;
  Py_ssize_t n = 0;
  std::vector<double> v;
  double operator[](Py_ssize_t i) const { return p[i]; }
};

// A sequence's items, borrowed from the list or tuple f.
PyObject *item(PyObject *f, Py_ssize_t i) {
  return PyList_Check(f) ? PyList_GetItem(f, i) : PyTuple_GetItem(f, i);
}
Py_ssize_t size(PyObject *f) { return PyList_Check(f) ? PyList_Size(f) : PyTuple_Size(f); }

// Numbers from o: an array of a numeric type (no copy for contiguous
// float64), a number, or a sequence of numbers and None.
Got readd(PyObject *o, Dv &d) {
  Ai a;
  if (arrif(o, a)) {
    if (numkind(a)) {
      d.n = a.n;
      if (a.big > 1) return MATRIX;
      if (a.kind == 'f' && a.size == 8 && (a.step == 8 || a.n <= 1) &&
          reinterpret_cast<std::uintptr_t>(a.data) % alignof(double) == 0) {
        d.p = reinterpret_cast<const double *>(a.data);
        return NUM;
      }
      d.v.resize(static_cast<std::size_t>(a.n));
      for (Py_ssize_t i = 0; i < a.n; ++i) d.v[static_cast<std::size_t>(i)] = aelt(a, i);
      d.p = d.v.data();
      return NUM;
    }
    // object arrays and other float sizes are read item by item below
    if (a.kind != 'O' && a.kind != 'f') return NOTNUM;
    if (a.big > 1) {
      d.n = a.n;
      return MATRIX;
    }
  }
  if (PyUnicode_Check(o) || PyBytes_Check(o) || PyByteArray_Check(o) || o == Py_None) {
    return NOTNUM;
  }
  if (!PySequence_Check(o)) {
    // a single number: Python has no vector of one
    d.v.resize(1);
    if (!PyNumber_Check(o) || !elt(o, d.v[0])) return NOTNUM;
    d.p = d.v.data();
    d.n = 1;
    return NUM;
  }
  PyObject *raw = PySequence_Fast(o, "not a sequence");
  if (raw == nullptr) {
    PyErr_Clear();
    return NOTNUM;
  }
  Obj f(raw);
  d.n = size(raw);
  d.v.resize(static_cast<std::size_t>(d.n));
  for (Py_ssize_t i = 0; i < d.n; ++i) {
    if (!elt(item(raw, i), d.v[static_cast<std::size_t>(i)])) return NOTNUM;
  }
  d.p = d.v.data();
  return NUM;
}

// TRUE, FALSE and NA (-1) read from o: a bool array, a bool, or a sequence
// of bools and None.
Got readb(PyObject *o, std::vector<int> &b) {
  Ai a;
  if (arrif(o, a)) {
    if (a.kind == 'b' && a.size == 1) {
      b.assign(static_cast<std::size_t>(a.n), 0);
      if (a.big > 1) return MATRIX;
      for (Py_ssize_t i = 0; i < a.n; ++i) b[static_cast<std::size_t>(i)] = a.data[i * a.step] != 0;
      return NUM;
    }
    if (a.kind != 'O') return NOTNUM;
    if (a.big > 1) {
      b.assign(static_cast<std::size_t>(a.n), 0);
      return MATRIX;
    }
  }
  if (isbool(o)) {
    const int t = PyObject_IsTrue(o);
    if (t < 0) throw Raised{};
    b.assign(1, t);
    return NUM;
  }
  if (PyUnicode_Check(o) || PyBytes_Check(o) || !PySequence_Check(o)) return NOTNUM;
  PyObject *raw = PySequence_Fast(o, "not a sequence");
  if (raw == nullptr) {
    PyErr_Clear();
    return NOTNUM;
  }
  Obj f(raw);
  const Py_ssize_t n = size(raw);
  b.assign(static_cast<std::size_t>(n), 0);
  for (Py_ssize_t i = 0; i < n; ++i) {
    PyObject *e = item(raw, i);
    if (e == Py_None) {
      b[static_cast<std::size_t>(i)] = -1;
      continue;
    }
    if (!isbool(e)) return NOTNUM;
    const int t = PyObject_IsTrue(e);
    if (t < 0) throw Raised{};
    b[static_cast<std::size_t>(i)] = t;
  }
  return NUM;
}

[[noreturn]] void notvec(const char *nm) {
  fail(false, "`%s` must be a vector, not a matrix with more than one row and column.", nm);
}

// .num(): finite numbers, of length len (unless len < 0), all at least lo
// (above lo when open).
void num(PyObject *o, const char *nm, Dv &d, Py_ssize_t len = -1,
         double lo = -HUGE_VAL, bool open = false) {
  const Got g = readd(o, d);
  if (g == NOTNUM || d.n < 1) {
    fail(g == NOTNUM, "`%s` must be finite numbers with no missing values.", nm);
  }
  if (g == MATRIX) notvec(nm);
  for (Py_ssize_t i = 0; i < d.n; ++i) {
    if (!std::isfinite(d[i])) fail(false, "`%s` must be finite numbers with no missing values.", nm);
  }
  if (len >= 0 && d.n != len) {
    fail(false, "`%s` must have length %zd, not %zd.", nm, len, d.n);
  }
  for (Py_ssize_t i = 0; i < d.n; ++i) {
    if (open ? d[i] <= lo : d[i] < lo) {
      fail(false, "`%s` must be %s%s.", nm, open ? "greater than " : "at least ", rfmt(lo).c_str());
    }
  }
}

// A per-group argument (src/bind.cpp, gvec() then .num()): not numbers at
// all is refused with the forms a group argument takes.
void gnum(PyObject *o, const char *nm, Dv &d, Py_ssize_t len = -1, double lo = -HUGE_VAL,
          bool open = false) {
  if (readd(o, d) == NOTNUM) {
    fail(true, "`%s` must be finite numbers, one per group: a numeric vector, a 1-d table, a "
               "list of single numbers or a one-row data frame.", nm);
  }
  d = Dv();
  num(o, nm, d, len, lo, open);
}

// x[1] when o is one number, NaN otherwise; type tells whether o was not
// numbers at all.
double val1(PyObject *o, bool &type) {
  Dv d;
  const Got g = readd(o, d);
  type = g == NOTNUM;
  return g == NUM && d.n == 1 ? d[0] : NaN;
}

// .one(): one number strictly between lo and hi.
double one(PyObject *o, const char *nm, double lo, double hi) {
  bool type;
  const double v = val1(o, type);
  if (std::isnan(v) || v <= lo || v >= hi) {
    fail(type, "`%s` must be one number strictly between %s and %s.", nm, rfmt(lo).c_str(),
         rfmt(hi).c_str());
  }
  return v;
}

// .p01(): one number from 0 to 1.
double p01(PyObject *o, const char *nm) {
  bool type;
  const double v = val1(o, type);
  if (std::isnan(v) || v < 0 || v > 1) fail(type, "`%s` must be one number from 0 to 1.", nm);
  return v;
}

// .whole(): one whole number from lo up to INT_MAX.
int whole(PyObject *o, const char *nm, int lo) {
  bool type;
  const double v = val1(o, type);
  if (std::isnan(v) || v != std::floor(v) || v < lo || v > INT_MAX) {
    fail(type, "`%s` must be one whole number from %d.", nm, lo);
  }
  return static_cast<int>(v);
}

// .lgl(): TRUE or FALSE, n of them (any number when n < 0), no missing
// value. Returns how many.
Py_ssize_t lgl(PyObject *o, const char *nm, Py_ssize_t n, std::vector<int> &b) {
  const Got g = readb(o, b);
  const Py_ssize_t k = static_cast<Py_ssize_t>(b.size());
  bool na = false;
  for (int v : b) na = na || v < 0;
  if (g == NOTNUM || (n >= 0 && k != n) || na) {
    fail(g == NOTNUM, "`%s` must be TRUE or FALSE, one per unit, with no missing values.", nm);
  }
  if (g == MATRIX) notvec(nm);
  return k;
}

// Building results.

PyObject *pyd(double x) {
  if (isna(x)) Py_RETURN_NONE;
  return PyFloat_FromDouble(x);
}

Obj flist(const double *x, Py_ssize_t n) {
  Obj l = own(PyList_New(n));
  for (Py_ssize_t i = 0; i < n; ++i) PyList_SetItem(l.get(), i, own(pyd(x[i])).release());
  return l;
}

Obj flist(const std::vector<double> &x) {
  return flist(x.data(), static_cast<Py_ssize_t>(x.size()));
}

Obj ilist(const std::vector<int> &x) {
  Obj l = own(PyList_New(static_cast<Py_ssize_t>(x.size())));
  for (std::size_t i = 0; i < x.size(); ++i) {
    PyList_SetItem(l.get(), static_cast<Py_ssize_t>(i), own(PyLong_FromLong(x[i])).release());
  }
  return l;
}

void put(const Obj &d, const char *k, Obj v) {
  if (PyDict_SetItemString(d.get(), k, v.get()) < 0) throw Raised{};
}

Obj dict() { return own(PyDict_New()); }

// A named vector: {name: value}.
Obj nvec(int n, const char *const *nm, const double *v) {
  Obj d = dict();
  for (int i = 0; i < n; ++i) put(d, nm[i], own(pyd(v[i])));
  return d;
}

// A column of bools, True at the reference group i (from 1) of n.
Obj refcol(int n, int i) {
  Obj l = own(PyList_New(n));
  for (int k = 0; k < n; ++k) PyList_SetItem(l.get(), k, own(PyBool_FromLong(k == i - 1)).release());
  return l;
}

// .grp(): the group names given, or "1", "2", ...
Obj grp(PyObject *names, Py_ssize_t n) {
  if (names != Py_None) return own(PySequence_List(names));
  Obj l = own(PyList_New(n));
  for (Py_ssize_t i = 0; i < n; ++i) {
    PyList_SetItem(l.get(), i, own(PyUnicode_FromFormat("%zd", i + 1)).release());
  }
  return l;
}

// .ref(): the group named or numbered (from 1) by r. A number is a
// position; None and bools are neither; anything else is a name.
int ref(PyObject *r, const Obj &nms, int n) {
  const bool pos = r == Py_None || isbool(r) || PyNumber_Check(r);
  if (!pos) {
    for (Py_ssize_t i = 0; i < n; ++i) {
      const int eq = PyObject_RichCompareBool(r, PyList_GetItem(nms.get(), i), Py_EQ);
      if (eq < 0) throw Raised{};
      if (eq) return static_cast<int>(i + 1);
    }
    Obj s = own(PyObject_Str(r));
    const char *c = PyUnicode_AsUTF8AndSize(s.get(), nullptr);
    if (c == nullptr) throw Raised{};
    fail(false, "`ref` names no group: %s.", c);
  }
  bool type;
  const double v = val1(r, type);
  if (std::isnan(v) || v != std::floor(v) || v < 1 || v > n) {
    fail(type, "`ref` must be one group name or a position from 1 to %d.", n);
  }
  return static_cast<int>(v);
}

// TRUE when some p[i] > 0 (with keep[i] when keep is given).
bool pos(const Dv &p, const std::vector<int> *keep = nullptr) {
  for (Py_ssize_t i = 0; i < p.n; ++i) {
    if ((keep == nullptr || (*keep)[static_cast<std::size_t>(i)]) && p[i] > 0) return true;
  }
  return false;
}

// rep(1, n)
void ones(Dv &d, Py_ssize_t n) {
  d.v.assign(static_cast<std::size_t>(n), 1.0);
  d.p = d.v.data();
  d.n = n;
}

// The entry points: a[] holds the positional arguments.

// Disparity.

PyObject *vodtst(PyObject *const *a) {
  who = "vodtst";
  Dv d, ni, sh;
  gnum(a[0], "day", d, -1, 0.0);
  const int n = ilen(static_cast<double>(d.n), "day");
  gnum(a[1], "night", ni, n, 0.0, true);
  const bool has = a[3] != Py_None;
  if (has) gnum(a[3], "shift", sh, n, 0.0, true);
  Obj nms = grp(a[4], n);
  const int i = ref(a[2], nms, n);
  if (d[i - 1] <= 0) fail(false, "the reference group has no daylight stops.");
  const std::size_t un = static_cast<std::size_t>(n);
  std::vector<double> sd(un), sn(un), ra(un), adj(has ? un : 0);
  std::vector<int> sg(un);
  kchk(clab_vod(d.p, ni.p, n, i - 1, has ? sh.p : nullptr, sd.data(), sn.data(), ra.data(),
                sg.data(), has ? adj.data() : nullptr));
  Obj out = dict();
  put(out, "group", std::move(nms));
  put(out, "ref", refcol(n, i));
  put(out, "shday", flist(sd));
  put(out, "shnight", flist(sn));
  put(out, "ratio", flist(ra));
  put(out, "sign", ilist(sg));
  if (has) put(out, "adj", flist(adj));
  return out.release();
}

PyObject *orrisk(PyObject *const *a) {
  who = "orrisk";
  Dv o;
  num(a[0], "or", o, 1, 0.0, true);
  const double v = o[0];
  if (!std::isfinite(1 / v)) fail(false, "`or` is too small: its reciprocal is not finite.");
  const double q = one(a[1], "rate", 0, 1);
  const double s = one(a[2], "share", 0, 1);
  double r[4];
  kchk(clab_orrisk(v, q, s, r, r + 1));
  r[2] = r[0] / r[1];
  r[3] = v;
  static const char *const nm[] = {"exposed", "unexposed", "rr", "or"};
  return nvec(4, nm, r).release();
}

PyObject *outtst(PyObject *const *a) {
  who = "outtst";
  Dv q, t, w;
  num(a[0], "q", q, -1, 0.0);
  const int n = ilen(static_cast<double>(q.n), "q");
  for (int i = 0; i < n; ++i) {
    if (q[i] > 1) fail(false, "`q` must be at most 1.");
  }
  num(a[1], "t", t, 1);
  if (a[2] == Py_None) ones(w, n); else num(a[2], "w", w, n, 0.0);
  if (!pos(w)) fail(false, "`w` must have a positive total.");
  double r[2];
  kchk(clab_outtst(q.p, w.p, n, t[0], r, r + 1));
  static const char *const nm[] = {"stop", "hit"};
  return nvec(2, nm, r).release();
}

PyObject *dspbm(PyObject *const *a) {
  who = "dspbm";
  Dv p, c, f;
  gnum(a[0], "pop", p, -1, 0.0, true);
  const int n = ilen(static_cast<double>(p.n), "pop");
  gnum(a[1], "contact", c, n, 0.0, true);
  gnum(a[2], "force", f, n, 0.0, true);
  Obj nms = grp(a[4], n);
  const int i = ref(a[3], nms, n);
  const std::size_t un = static_cast<std::size_t>(n);
  std::vector<double> cd(un), fd(un), rdv(un);
  kchk(clab_dspbm(p.p, c.p, f.p, n, i - 1, cd.data(), fd.data(), rdv.data()));
  Obj out = dict();
  put(out, "group", std::move(nms));
  put(out, "ref", refcol(n, i));
  put(out, "cd", flist(cd));
  put(out, "fd", flist(fd));
  put(out, "rd", flist(rdv));
  return out.release();
}

// a: count (its values column by column, or None when it is not a table),
// its rows, its columns, pop, the group names.
PyObject *dyadrt(PyObject *const *a) {
  who = "dyadrt";
  Dv p, x;
  gnum(a[3], "pop", p, -1, 0.0, true);
  const int g = ilen(static_cast<double>(p.n), "pop");
  Obj nms = grp(a[4], g);
  const long nr = PyLong_AsLong(a[1]), nc = PyLong_AsLong(a[2]);
  if (PyErr_Occurred()) throw Raised{};
  if (a[0] == Py_None || nr != g || nc != g) fail(false, "`count` must be a %d by %d matrix.", g, g);
  num(a[0], "count", x, -1, 0.0);
  const Py_ssize_t m = x.n;
  if (!pos(x)) fail(false, "`count` holds no offences.");
  const std::size_t um = static_cast<std::size_t>(m);
  std::vector<int> off(um), vic(um);
  Obj so = own(PyList_New(m)), sv = own(PyList_New(m));
  for (Py_ssize_t k = 0; k < m; ++k) {
    off[static_cast<std::size_t>(k)] = static_cast<int>(k % g);
    vic[static_cast<std::size_t>(k)] = static_cast<int>(k / g);
    PyObject *o = PyList_GetItem(nms.get(), k % g), *v = PyList_GetItem(nms.get(), k / g);
    Py_INCREF(o);
    PyList_SetItem(so.get(), k, o);
    Py_INCREF(v);
    PyList_SetItem(sv.get(), k, v);
  }
  std::vector<double> ro(um), rp(um), no(um), ra(um);
  kchk(clab_dyadrt(x.p, off.data(), vic.data(), ilen(static_cast<double>(m), "count"), p.p, g,
                   ro.data(), rp.data(), no.data(), ra.data()));
  Obj out = dict();
  put(out, "off", std::move(so));
  put(out, "vic", std::move(sv));
  put(out, "count", flist(x.p, m));
  put(out, "rtoff", flist(ro));
  put(out, "rtpair", flist(rp));
  put(out, "nulloff", flist(no));
  put(out, "ratio", flist(ra));
  return out.release();
}

// Forensics.

PyObject *gtmp(PyObject *const *a) {
  who = "gtmp";
  Dv p;
  gnum(a[0], "p", p, -1, 0.0, true);
  const int n = ilen(static_cast<double>(p.n), "p");
  for (int i = 0; i < n; ++i) {
    if (p[i] > 1) fail(false, "`p` must be at most 1.");
  }
  std::vector<int> h0;
  const Got gh = readb(a[2], h0);
  const Py_ssize_t nh = static_cast<Py_ssize_t>(h0.size());
  bool na = false;
  for (int v : h0) na = na || v < 0;
  if (gh == NOTNUM || na || (nh != 1 && nh != n)) {
    fail(gh == NOTNUM, "`hom` must be TRUE or FALSE, once or once per locus.");
  }
  if (gh == MATRIX) notvec("hom");
  const std::size_t un = static_cast<std::size_t>(n);
  std::vector<int> h(un);
  for (int i = 0; i < n; ++i) h[static_cast<std::size_t>(i)] = h0[nh == 1 ? 0 : static_cast<std::size_t>(i)];
  std::vector<double> q(un, NaN);
  if (a[1] != Py_None) {
    Dv qd;
    const Got gq = readd(a[1], qd);
    if (gq == NOTNUM || qd.n != n) fail(gq == NOTNUM, "`q` must be numbers, one per locus.");
    if (gq == MATRIX) notvec("q");
    for (int i = 0; i < n; ++i) q[static_cast<std::size_t>(i)] = qd[i];
  }
  bool bad = false;
  for (std::size_t i = 0; i < un; ++i) {
    if (h[i]) q[i] = 0;
    const double pi = p[static_cast<Py_ssize_t>(i)];
    bad = bad || std::isnan(q[i]) || (!h[i] && (q[i] <= 0 || pi + q[i] > 1 + 1e-12));
  }
  if (bad) {
    fail(false, "`q` must be above 0 with `p + q` at most 1 at every heterozygous locus.");
  }
  bool type;
  const double th = val1(a[3], type);
  if (std::isnan(th) || th < 0 || th >= 1) {
    fail(type, "`theta` must be one number from 0 up to, not including, 1.");
  }
  std::vector<double> mp(un);
  double l10 = 0;
  kchk(clab_gtmp(p.p, q.data(), h.data(), n, th, mp.data(), &l10));
  Obj loc = dict();
  put(loc, "locus", grp(a[4], n));
  Obj hl = own(PyList_New(n));
  for (int i = 0; i < n; ++i) {
    PyList_SetItem(hl.get(), i, own(PyBool_FromLong(h[static_cast<std::size_t>(i)])).release());
  }
  put(loc, "hom", std::move(hl));
  put(loc, "mp", flist(mp));
  Obj out = dict();
  put(out, "locus", std::move(loc));
  put(out, "log10mp", own(PyFloat_FromDouble(l10)));
  put(out, "log10lr", own(PyFloat_FromDouble(-l10)));
  return out.release();
}

PyObject *cllr(PyObject *const *a) {
  who = "cllr";
  Dv s, d;
  num(a[0], "same", s, -1, 0.0, true);
  num(a[1], "diff", d, -1, 0.0, true);
  double r = 0;
  kchk(clab_cllr(s.p, ilen(static_cast<double>(s.n), "same"), d.p,
                 ilen(static_cast<double>(d.n), "diff"), &r));
  return PyFloat_FromDouble(r);
}

PyObject *lrpost(PyObject *const *a) {
  who = "lrpost";
  Dv l, p;
  num(a[0], "lr", l, -1, 0.0, true);
  num(a[1], "prior", p, -1, 0.0, true);
  bool bad = p.n != 1 && p.n != l.n;
  for (Py_ssize_t i = 0; i < p.n && !bad; ++i) bad = p[i] >= 1;
  if (bad) fail(false, "`prior` must be below 1, once or once per ratio.");
  std::vector<double> o(static_cast<std::size_t>(l.n));
  for (Py_ssize_t i = 0; i < l.n; ++i) {
    const double q = p[p.n == 1 ? 0 : i];
    o[static_cast<std::size_t>(i)] = 1 / (1 + (1 - q) / (q * l[i]));
  }
  return flist(o).release();
}

// A set of the items of a list; null (no error set) when one cannot be
// hashed.
Obj set(PyObject *l) {
  PyObject *s = PySet_New(l);
  if (s == nullptr) PyErr_Clear();
  return Obj(s);
}

// .mix(): one side of a mixture ratio.
double mix(const Dv &p, PyObject *names, const Obj &nset, PyObject *need, PyObject *x,
           const char *an, const char *bn) {
  bool bad = !PyList_Check(need);
  Obj ns;
  if (!bad) {
    const Py_ssize_t k = PyList_Size(need);
    for (Py_ssize_t i = 0; i < k && !bad; ++i) {
      PyObject *e = PyList_GetItem(need, i);
      bad = e == Py_None;
      if (!bad) {
        const int in = PySet_Contains(nset.get(), e);
        if (in < 0) PyErr_Clear();
        bad = in != 1;
      }
    }
    if (!bad) {
      ns = set(need);
      bad = PySet_Size(ns.get()) > 20;
    }
  }
  if (bad) fail(false, "`%s` must be names of `p`, at most 20.", an);
  bool type;
  const double v = val1(x, type);
  if (std::isnan(v) || v < 0 || v != std::floor(v) || v > 50) {
    fail(type, "`%s` must be a whole number from 0 to 50.", bn);
  }
  std::vector<int> in(static_cast<std::size_t>(p.n));
  for (Py_ssize_t i = 0; i < p.n; ++i) {
    const int c = PySet_Contains(ns.get(), PyList_GetItem(names, i));
    if (c < 0) throw Raised{};
    in[static_cast<std::size_t>(i)] = c;
  }
  double r = 0;
  kchk(clab_mixp(p.p, in.data(), static_cast<int>(p.n), static_cast<int>(v), &r));
  return r;
}

// a: p, its names, np, xp, nd, xd
PyObject *mixlr(PyObject *const *a) {
  who = "mixlr";
  Dv p;
  num(a[0], "p", p, -1, 0.0, true);
  const int n = ilen(static_cast<double>(p.n), "p");
  bool bad = a[1] == Py_None;
  Obj nms, nset;
  if (!bad) {
    nms = own(PySequence_List(a[1]));
    nset = set(nms.get());
    bad = nset.get() == nullptr || PySet_Size(nset.get()) != n;
  }
  if (bad || nsum(p.p, n) > 1 + 1e-12) {
    fail(false, "`p` must be uniquely named with a sum of at most 1.");
  }
  double r[3];
  r[0] = mix(p, nms.get(), nset, a[2], a[3], "np", "xp");
  r[1] = mix(p, nms.get(), nset, a[4], a[5], "nd", "xd");
  r[2] = r[0] / r[1];
  static const char *const nm[] = {"num", "den", "lr"};
  return nvec(3, nm, r).release();
}

PyObject *dbsrch(PyObject *const *a) {
  who = "dbsrch";
  Dv m, nn;
  num(a[0], "mp", m, 1, 0.0, true);
  const double mv = m[0];
  if (mv > 1) fail(false, "`mp` must be at most 1.");
  num(a[1], "n", nn, 1, 1.0);
  const double nv = nn[0];
  double r[3];
  const double e = nv * mv;
  r[0] = e < 1 ? e : 1;
  r[1] = -std::expm1(nv * std::log1p(-mv));
  const bool has = a[2] != Py_None;
  if (has) {
    Dv pp;
    num(a[2], "pop", pp, 1, nv);
    r[2] = 1 / (1 + rd((pp[0] - nv) * mv));
  }
  static const char *const nm[] = {"np", "any", "post"};
  return nvec(has ? 3 : 2, nm, r).release();
}

// .cal(): the ratios checked, and calibrated in the order given.
struct Cal {
  Dv s, d;
  std::vector<double> cs, cd;
  double min = 0;
};

void cal(PyObject *same, PyObject *diff, Cal &c) {
  num(same, "same", c.s, -1, 0.0, true);
  num(diff, "diff", c.d, -1, 0.0, true);
  const int ns = ilen(static_cast<double>(c.s.n), "same");
  const int nd = ilen(static_cast<double>(c.d.n), "diff");
  const int n = ilen(static_cast<double>(ns) + nd, "same");
  const std::size_t un = static_cast<std::size_t>(n);
  // lr <- c(same, diff); o <- order(lr), stable as R's radix order
  std::vector<double> lr(un);
  std::copy(c.s.p, c.s.p + ns, lr.begin());
  std::copy(c.d.p, c.d.p + nd, lr.begin() + ns);
  std::vector<int> o(un);
  for (int k = 0; k < n; ++k) o[static_cast<std::size_t>(k)] = k;
  std::stable_sort(o.begin(), o.end(), [&lr](int i, int j) {
    return lr[static_cast<std::size_t>(i)] < lr[static_cast<std::size_t>(j)];
  });
  std::vector<int> is(un);
  std::vector<double> ls(un), post(un), work(2 * un);
  for (std::size_t k = 0; k < un; ++k) {
    is[k] = o[k] < ns;
    ls[k] = lr[static_cast<std::size_t>(o[k])];
  }
  kchk(clab_lrcal(ls.data(), is.data(), n, post.data(), work.data(), &c.min));
  c.cs.resize(static_cast<std::size_t>(ns));
  c.cd.resize(static_cast<std::size_t>(nd));
  for (std::size_t k = 0; k < un; ++k) {
    const double x = post[k] / (1 - post[k]) * nd / ns;
    if (o[k] < ns) c.cs[static_cast<std::size_t>(o[k])] = x;
    else c.cd[static_cast<std::size_t>(o[k] - ns)] = x;
  }
}

PyObject *cllrmin(PyObject *const *a) {
  who = "cllrmin";
  Cal c;
  cal(a[0], a[1], c);
  double r[3];
  kchk(clab_cllr(c.s.p, static_cast<int>(c.s.n), c.d.p, static_cast<int>(c.d.n), r));
  r[1] = c.min;
  r[2] = r[0] - c.min;
  static const char *const nm[] = {"cllr", "min", "cal"};
  return nvec(3, nm, r).release();
}

PyObject *lrcal(PyObject *const *a) {
  who = "lrcal";
  Cal c;
  cal(a[0], a[1], c);
  Obj out = dict();
  put(out, "same", flist(c.cs));
  put(out, "diff", flist(c.cd));
  return out.release();
}

// one curve of tipplt(): sorted log10 ratios and the share at or above
Obj curve(const Dv &x) {
  const std::size_t n = static_cast<std::size_t>(x.n);
  std::vector<double> l(n), s(n);
  for (std::size_t i = 0; i < n; ++i) l[i] = std::log10(x[static_cast<Py_ssize_t>(i)]);
  std::sort(l.begin(), l.end());
  for (std::size_t i = 0; i < n; ++i) s[i] = static_cast<double>(n - i) / static_cast<double>(n);
  Obj out = dict();
  put(out, "log10lr", flist(l));
  put(out, "share", flist(s));
  return out;
}

PyObject *tipplt(PyObject *const *a) {
  who = "tipplt";
  Dv s, d;
  num(a[0], "same", s, -1, 0.0, true);
  num(a[1], "diff", d, -1, 0.0, true);
  ilen(static_cast<double>(s.n), "same");
  ilen(static_cast<double>(d.n), "diff");
  Obj out = dict();
  put(out, "same", curve(s));
  put(out, "diff", curve(d));
  return out.release();
}

// log(1 + e^z), without overflow for large z
double softplus(double z) { return z > 0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z)); }

// One empirical cross-entropy, as src/bind.cpp computes it.
double ece(const double *s, Py_ssize_t ns, const double *d, Py_ssize_t nd, double lo,
           double odds, double pr, double *t) {
  const double lg = lo * std::log(10.0);
  const bool fin = std::isfinite(odds) && odds > 0;
  for (Py_ssize_t k = 0; k < ns; ++k) {
    const double x = fin && std::isfinite(s[k]) && s[k] > 0 ? s[k] * odds : 0;
    if (std::isfinite(x) && x > 0) {
      t[k] = x < 1 ? std::log1p(x) - std::log(x) : std::log1p(1 / x);
    } else {
      t[k] = softplus(-(std::log(s[k]) + lg));
    }
  }
  const double a = nmean(t, ns);
  for (Py_ssize_t k = 0; k < nd; ++k) {
    const double x = fin && std::isfinite(d[k]) && d[k] > 0 ? d[k] * odds : 0;
    t[k] = std::isfinite(x) && x > 0 ? std::log1p(x) : softplus(std::log(d[k]) + lg);
  }
  const double b = nmean(t, nd);
  return (rd(pr * a) + rd((1 - pr) * b)) / std::log(2.0);
}

PyObject *eceplt(PyObject *const *a) {
  who = "eceplt";
  Cal c;
  cal(a[0], a[1], c);
  Dv l;
  num(a[2], "lo", l);
  const Py_ssize_t nl = l.n, ns = c.s.n, nd = c.d.n;
  const std::size_t ul = static_cast<std::size_t>(nl);
  std::vector<double> e(ul), m(ul), r(ul), t(static_cast<std::size_t>(ns > nd ? ns : nd));
  const double one = 1.0;
  for (Py_ssize_t i = 0; i < nl; ++i) {
    const double odds = rpow(10.0, l[i]);
    const double pr = std::isinf(odds) ? 1.0 : odds / (1 + odds);
    const std::size_t ui = static_cast<std::size_t>(i);
    e[ui] = ece(c.s.p, ns, c.d.p, nd, l[i], odds, pr, t.data());
    m[ui] = ece(c.cs.data(), ns, c.cd.data(), nd, l[i], odds, pr, t.data());
    r[ui] = ece(&one, 1, &one, 1, l[i], odds, pr, t.data());
  }
  Obj out = dict();
  put(out, "lo", flist(l.p, nl));
  put(out, "ece", flist(e));
  put(out, "min", flist(m));
  put(out, "ref", flist(r));
  return out.release();
}

// Whole numbers from a list Python built (allele codes).
std::vector<int> ints(PyObject *l) {
  const Py_ssize_t n = PyList_Size(l);
  if (n < 0) throw Raised{};
  std::vector<int> v(static_cast<std::size_t>(n));
  for (Py_ssize_t i = 0; i < n; ++i) {
    v[static_cast<std::size_t>(i)] = static_cast<int>(PyLong_AsLong(PyList_GetItem(l, i)));
  }
  if (PyErr_Occurred()) throw Raised{};
  return v;
}

// .gf(): the frequencies of the alleles of a table of nr rows whose codes
// are g (two per row), given row by row; None when freq was not a table of
// that shape. gn names the table in the message.
void gf(PyObject *freq, const std::vector<int> &g, Py_ssize_t nr, const char *gn, Dv &f) {
  Got got = NOTNUM;
  if (freq != Py_None) got = readd(freq, f);
  bool bad = got != NUM || f.n != 2 * nr;
  for (Py_ssize_t i = 0; !bad && i < f.n; ++i) bad = std::isnan(f[i]) || f[i] <= 0 || f[i] > 1;
  if (bad) {
    fail(got == NOTNUM && freq != Py_None,
         "`freq` must hold the frequencies of the alleles of `%s`, same shape, above 0 and at "
         "most 1.",
         gn);
  }
  for (Py_ssize_t i = 0; i < nr; ++i) {
    const std::size_t k = static_cast<std::size_t>(2 * i);
    if (g[k] == g[k + 1] && f[2 * i] != f[2 * i + 1]) {
      fail(false, "`freq` gives two frequencies for one allele of a homozygote.");
    }
  }
}

// .kin(): the identity-by-descent probabilities of a relationship.
void kin(PyObject *rel, double *k) {
  static const char *const nm[] = {"parent", "sibling", "half", "cousin"};
  static const double ks[4][3] = {{0, 1, 0}, {0.25, 0.5, 0.25}, {0.5, 0.5, 0}, {0.75, 0.25, 0}};
  if (PyUnicode_Check(rel)) {
    const char *s = PyUnicode_AsUTF8AndSize(rel, nullptr);
    if (s == nullptr) throw Raised{};
    for (int i = 0; i < 4; ++i) {
      if (std::strcmp(s, nm[i]) == 0) {
        std::memcpy(k, ks[i], sizeof ks[i]);
        return;
      }
    }
  } else {
    Dv v;
    if (readd(rel, v) == NUM && v.n == 3) {
      bool ok = true;
      for (int i = 0; i < 3; ++i) ok = ok && !std::isnan(v[i]) && v[i] >= 0;
      if (ok && std::fabs(nsum(v.p, 3) - 1) < 1e-9) {
        std::memcpy(k, v.p, 3 * sizeof(double));
        return;
      }
    }
  }
  fail(false, "`rel` must be one of parent, sibling, half, cousin, or three probabilities that "
              "sum to 1.");
}

// .idx(): the per-locus indices with their loci, their log10 total and the
// posterior at even prior odds.
PyObject *idx(PyObject *loci, const std::vector<double> &ix) {
  const std::size_t nr = ix.size();
  std::vector<double> lg(nr);
  for (std::size_t i = 0; i < nr; ++i) lg[i] = std::log10(ix[i]);
  const double l10 = nsum(lg.data(), static_cast<Py_ssize_t>(nr));
  Obj loc = dict();
  put(loc, "locus", grp(loci, static_cast<Py_ssize_t>(nr)));
  put(loc, "idx", flist(ix));
  Obj out = dict();
  put(out, "locus", std::move(loc));
  put(out, "log10idx", own(PyFloat_FromDouble(l10)));
  put(out, "post", own(PyFloat_FromDouble(1 / (1 + rpow(10.0, -l10)))));
  return out.release();
}

// a: codes of g1 and g2 (two per row), freq row by row (or None), rel,
// the loci (or None), the name of g1 for messages
PyObject *kinidx(PyObject *const *a) {
  who = "kinidx";
  const std::vector<int> g1 = ints(a[0]), g2 = ints(a[1]);
  const Py_ssize_t nr = static_cast<Py_ssize_t>(g1.size() / 2);
  Dv f;
  gf(a[2], g1, nr, "g1", f);
  double k[3];
  kin(a[3], k);
  std::vector<double> ix(static_cast<std::size_t>(nr));
  kchk(clab_kinidx(g1.data(), g2.data(), f.p, static_cast<int>(nr), k, ix.data()));
  return idx(a[4], ix);
}

// a: codes of child, father and mother (or None), freq, the loci
PyObject *patidx(PyObject *const *a) {
  who = "patidx";
  const std::vector<int> c = ints(a[0]), fa = ints(a[1]);
  const Py_ssize_t nr = static_cast<Py_ssize_t>(c.size() / 2);
  Dv f;
  std::vector<double> ix(static_cast<std::size_t>(nr));
  if (a[2] == Py_None) {
    // kinidx(child, father, freq, "parent"), its messages naming g1
    gf(a[3], c, nr, "g1", f);
    const double k[3] = {0, 1, 0};
    kchk(clab_kinidx(c.data(), fa.data(), f.p, static_cast<int>(nr), k, ix.data()));
  } else {
    const std::vector<int> m = ints(a[2]);
    gf(a[3], c, nr, "child", f);
    kchk(clab_patidx(c.data(), m.data(), fa.data(), f.p, static_cast<int>(nr), ix.data()));
  }
  return idx(a[4], ix);
}

// Policing.

PyObject *clrsmp(PyObject *const *a) {
  who = "clrsmp";
  Dv ca, ra, cb, rb;
  gnum(a[0], "clra", ca, -1, 0.0);
  const int n = ilen(static_cast<double>(ca.n), "clra");
  gnum(a[1], "reca", ra, n, 0.0, true);
  gnum(a[2], "clrb", cb, n, 0.0);
  gnum(a[3], "recb", rb, n, 0.0, true);
  const std::size_t un = static_cast<std::size_t>(n);
  std::vector<double> c[6];
  for (auto &v : c) v.resize(un);
  double tot[8];
  int cnt[3];
  kchk(clab_clrsmp(ca.p, ra.p, cb.p, rb.p, n, c[0].data(), c[1].data(), c[2].data(),
                   c[3].data(), c[4].data(), c[5].data(), tot, cnt));
  static const char *const pnm[] = {"ratea", "rateb", "sharea", "shareb", "rateeff", "mixeff"};
  Obj parts = dict();
  put(parts, "part", grp(a[4], n));
  for (int k = 0; k < 6; ++k) put(parts, pnm[k], flist(c[k]));
  for (int k = 0; k < 3; ++k) tot[5 + k] = cnt[k];
  static const char *const tnm[] = {"ratea", "rateb", "change", "rateeff",
                                    "mixeff", "up",    "down",   "simpson"};
  Obj out = dict();
  put(out, "parts", std::move(parts));
  put(out, "total", nvec(8, tnm, tot));
  return out.release();
}

PyObject *expbnd(PyObject *const *a) {
  who = "expbnd";
  Dv y, m, g;
  gnum(a[0], "y", y, -1, 0.0, true);
  const int n = ilen(static_cast<double>(y.n), "y");
  gnum(a[1], "m", m, n, 0.0, true);
  num(a[2], "gamma", g, 1, 1.0);
  Obj nms = grp(a[4], n);
  const int i = ref(a[3], nms, n);
  const std::size_t un = static_cast<std::size_t>(n);
  std::vector<double> c[6];
  for (auto &v : c) v.resize(un);
  std::vector<int> sg(un);
  kchk(clab_expbnd(y.p, m.p, n, i - 1, g[0], c[0].data(), c[1].data(), c[2].data(),
                   c[3].data(), c[4].data(), c[5].data(), sg.data()));
  static const char *const nm[] = {"rate", "lo", "hi", "ratio", "rlo", "rhi"};
  Obj out = dict();
  put(out, "group", std::move(nms));
  put(out, "ref", refcol(n, i));
  for (int k = 0; k < 6; ++k) put(out, nm[k], flist(c[k]));
  put(out, "sign", ilist(sg));
  return out.release();
}

PyObject *misrd(PyObject *const *a) {
  who = "misrd";
  const double r1 = p01(a[0], "r1");
  const double r0 = p01(a[1], "r0");
  const double p = one(a[2], "prev", 0, 1);
  const double s = p01(a[3], "se");
  const double c = p01(a[4], "sp");
  const double d1 = rd(p * s) + rd((1 - p) * (1 - c));
  if (d1 <= 0 || d1 >= 1) fail(false, "`se` and `sp` leave a recorded group empty.");
  double r[5];
  kchk(clab_misrd(r1, r0, p, s, c, r));
  static const char *const nm[] = {"true", "measured", "factor", "d1", "d0"};
  return nvec(5, nm, r).release();
}

PyObject *rtmhs(PyObject *const *a) {
  who = "rtmhs";
  Dv x1, x2, c, w;
  num(a[0], "x1", x1);
  const int n = ilen(static_cast<double>(x1.n), "x1");
  num(a[1], "x2", x2, n);
  num(a[2], "cut", c, 1);
  if (a[3] == Py_None) ones(w, n); else num(a[3], "w", w, n, 0.0);
  bool any = false;
  for (int i = 0; i < n && !any; ++i) any = x1[i] > c[0] && w[i] > 0;
  if (!any) fail(false, "no place with weight exceeds `cut`.");
  double r[7];
  kchk(clab_rtmhs(x1.p, x2.p, w.p, n, c[0], r));
  static const char *const nm[] = {"nsel", "change", "mirror", "sym", "excess", "low", "lowsym"};
  return nvec(7, nm, r).release();
}

// Corrections.

PyObject *incap(PyObject *const *a) {
  who = "incap";
  Dv l, q, s, sh;
  gnum(a[0], "lam", l, -1, 0.0, true);
  const int n = ilen(static_cast<double>(l.n), "lam");
  gnum(a[1], "q", q, -1, 0.0, true);
  for (Py_ssize_t i = 0; i < q.n; ++i) {
    if (q[i] > 1) fail(false, "`q` must be at most 1.");
  }
  gnum(a[2], "s", s, -1, 0.0);
  if (q.n != 1 && q.n != n) fail(false, "`q` must have length 1 or %d.", n);
  if (s.n != 1 && s.n != n) fail(false, "`s` must have length 1 or %d.", n);
  const bool has = a[3] != Py_None;
  if (has) {
    gnum(a[3], "share", sh, n, 0.0);
    if (std::fabs(nsum(sh.p, n) - 1) > 1e-9) fail(false, "`share` must sum to 1.");
  }
  const std::size_t nr = static_cast<std::size_t>(has ? n + 1 : n);
  std::vector<double> c[6];
  for (auto &v : c) v.resize(nr);
  for (int i = 0; i < n; ++i) {
    const std::size_t u = static_cast<std::size_t>(i);
    c[0][u] = l[i];
    c[1][u] = q[q.n == 1 ? 0 : i];
    c[2][u] = s[s.n == 1 ? 0 : i];
  }
  double agg[3];
  kchk(clab_incap(c[0].data(), c[1].data(), c[2].data(), n, has ? sh.p : nullptr, c[3].data(),
                  c[4].data(), c[5].data(), agg));
  Obj g = grp(a[4], n);
  if (has) {
    const std::size_t u = static_cast<std::size_t>(n);
    c[0][u] = agg[0];
    c[1][u] = NA;
    c[2][u] = NA;
    c[3][u] = agg[1];
    c[4][u] = agg[2];
    c[5][u] = NA;
    Obj all = own(PyUnicode_FromString("all"));
    if (PyList_Append(g.get(), all.get()) < 0) throw Raised{};
  }
  static const char *const nm[] = {"lam", "q", "s", "rate", "prev", "marg"};
  Obj out = dict();
  put(out, "group", std::move(g));
  for (int k = 0; k < 6; ++k) put(out, nm[k], flist(c[k]));
  return out.release();
}

PyObject *incapr(PyObject *const *a) {
  who = "incapr";
  Dv l;
  num(a[0], "lam", l, -1, 0.0);
  for (Py_ssize_t i = 1; i < l.n; ++i) {
    if (l[i] - l[i - 1] > 0) fail(false, "`lam` must never increase.");
  }
  const int sv = whole(a[2], "s", 1);
  const int tv = whole(a[1], "t0", 0);
  if (static_cast<double>(tv) + sv > static_cast<double>(l.n) - 1) {
    fail(false, "`lam` must have length at least t0 + s + 1.");
  }
  const double rv = p01(a[3], "r");
  double o[5];
  kchk(clab_incapr(l.p, ilen(static_cast<double>(l.n), "lam"), tv, sv, rv, o));
  static const char *const nm[] = {"prev", "net", "upper", "lower", "later"};
  return nvec(5, nm, o).release();
}

PyObject *sellab(PyObject *const *a) {
  who = "sellab";
  std::vector<int> rel, rule;
  const Py_ssize_t n = lgl(a[1], "rel", -1, rel);
  lgl(a[2], "rule", n, rule);
  Dv y;
  const Got gy = readd(a[0], y);
  if (gy == NOTNUM || y.n != n) fail(gy == NOTNUM, "`y` must be numbers, one per defendant.");
  if (gy == MATRIX) notvec("y");
  std::vector<double> yv(y.p, y.p + n);
  for (Py_ssize_t i = 0; i < n; ++i) {
    const std::size_t u = static_cast<std::size_t>(i);
    if (rel[u] && !(yv[u] == 0 || yv[u] == 1)) {
      fail(false, "`y` must be 0 or 1 for every released defendant.");
    }
  }
  for (std::size_t i = 0; i < yv.size(); ++i) {
    if (!rel[i]) yv[i] = 0;
  }
  Dv w;
  if (a[3] == Py_None) ones(w, n); else num(a[3], "w", w, n, 0.0);
  if (!pos(w, &rel)) fail(false, "the judge released nobody.");
  if (!pos(w, &rule)) fail(false, "the rule releases nobody.");
  double r[6];
  kchk(clab_sellab(yv.data(), rel.data(), rule.data(), w.p, ilen(static_cast<double>(n), "rel"),
                   r));
  static const char *const nm[] = {"obs", "unobs", "lower", "upper", "width", "naive"};
  return nvec(6, nm, r).release();
}

PyObject *losbnd(PyObject *const *a) {
  who = "losbnd";
  Dv d, o;
  num(a[0], "done", d, -1, 0.0);
  num(a[1], "open", o, -1, 0.0);
  double r[4];
  kchk(clab_losbnd(d.p, ilen(static_cast<double>(d.n), "done"), o.p,
                   ilen(static_cast<double>(o.n), "open"), r));
  static const char *const nm[] = {"mean", "age", "lower", "bias"};
  return nvec(4, nm, r).release();
}

int eq(PyObject *x, PyObject *y) {
  const int e = PyObject_RichCompareBool(x, y, Py_EQ);
  if (e < 0) throw Raised{};
  return e;
}

// a: y, the labels of z as a list of str (None where missing), the label
// of treat (None for the default; any other object never matches), w, and
// the kind of z: "l" (bools), "n" (numbers) or "s" (anything else)
PyObject *sntbnd(PyObject *const *a) {
  who = "sntbnd";
  Dv y;
  num(a[0], "y", y, -1, 0.0);
  const int n = ilen(static_cast<double>(y.n), "y");
  for (int i = 0; i < n; ++i) {
    if (!(y[i] == 0 || y[i] == 1)) fail(false, "`y` must be 0 or 1.");
  }
  PyObject *z = a[1];
  bool bad = PyList_Size(z) != n;
  for (Py_ssize_t i = 0; i < n && !bad; ++i) bad = PyList_GetItem(z, i) == Py_None;
  if (bad) fail(false, "`z` must have length %d and no missing values.", n);
  // lv <- sort(unique(z)): exactly two labels, in code-point order
  PyObject *lv[2] = {PyList_GetItem(z, 0), nullptr};
  bool more = false;
  for (Py_ssize_t i = 1; i < n && !more; ++i) {
    PyObject *e = PyList_GetItem(z, i);
    if (eq(e, lv[0])) continue;
    if (lv[1] == nullptr) lv[1] = e;
    else more = !eq(e, lv[1]);
  }
  if (lv[1] == nullptr || more) fail(false, "`z` must take exactly two values.");
  const int c = PyUnicode_Compare(lv[1], lv[0]);
  if (c == -1 && PyErr_Occurred()) throw Raised{};
  if (c < 0) std::swap(lv[0], lv[1]);
  // treat: one of lv, never chosen by sorting; by default TRUE for a
  // logical z and 1 for a numeric z of 0 and 1, otherwise required
  PyObject *tr = lv[1];
  if (a[2] == Py_None) {
    const char *k = PyUnicode_AsUTF8AndSize(a[4], nullptr);
    const char *s0 = PyUnicode_AsUTF8AndSize(lv[0], nullptr);
    const char *s1 = PyUnicode_AsUTF8AndSize(lv[1], nullptr);
    if (k == nullptr || s0 == nullptr || s1 == nullptr) throw Raised{};
    const bool b = k[0] == 'n' && std::strcmp(s0, "0") == 0 && std::strcmp(s1, "1") == 0;
    if (k[0] != 'l' && !b) {
      fail(false,
           "`treat` must name the treatment, %s or %s: it is chosen by default only for a "
           "logical `z` (TRUE) or a `z` of 0 and 1 (1).",
           s0, s1);
    }
  } else {
    tr = a[2];
    if (!eq(tr, lv[0]) && !eq(tr, lv[1])) fail(false, "`treat` must be one of the values of `z`.");
    tr = eq(tr, lv[0]) ? lv[0] : lv[1];
  }
  Dv w;
  if (a[3] == Py_None) ones(w, n); else num(a[3], "w", w, n, 0.0);
  if (!pos(w)) fail(false, "`w` must have a positive total.");
  std::vector<int> zi(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) zi[static_cast<std::size_t>(i)] = eq(PyList_GetItem(z, i), tr);
  double r[11];
  kchk(clab_sntbnd(y.p, zi.data(), w.p, n, r));
  Obj term = own(PyList_New(3));
  PyObject *other = tr == lv[0] ? lv[1] : lv[0];
  Py_INCREF(other);
  PyList_SetItem(term.get(), 0, other);
  Py_INCREF(tr);
  PyList_SetItem(term.get(), 1, tr);
  PyList_SetItem(term.get(), 2, own(PyUnicode_FromString("effect")).release());
  const double col[4][3] = {{r[2], r[3], NA},
                            {r[4], r[6], r[8]},
                            {r[5], r[7], r[9]},
                            {r[0] / r[2], r[1] / r[3], r[10]}};
  static const char *const nm[] = {"share", "lower", "upper", "observed"};
  Obj out = dict();
  put(out, "term", std::move(term));
  for (int k = 0; k < 4; ++k) put(out, nm[k], flist(col[k], 3));
  return out.release();
}

PyObject *hzdsel(PyObject *const *a) {
  who = "hzdsel";
  const double s = one(a[0], "s", 0, 1);
  const double h = p01(a[1], "h");
  const double l = p01(a[2], "l");
  if (l >= h) fail(false, "`l` must be below `h`.");
  Dv sa, sb;
  gnum(a[3], "surh", sa, 2, 0.0, true);
  gnum(a[4], "surl", sb, 2, 0.0, true);
  for (int i = 0; i < 2; ++i) {
    if (sa[i] > 1) fail(false, "`surh` must be at most 1.");
  }
  for (int i = 0; i < 2; ++i) {
    if (sb[i] > 1) fail(false, "`surl` must be at most 1.");
  }
  double r[5];
  kchk(clab_hzdsel(s, h, l, sa.p, sb.p, r));
  static const char *const nm[] = {"share0", "share1", "haz0", "haz1", "ratio"};
  return nvec(5, nm, r).release();
}

// The self-test (src/selftest.cpp): checks run, checks failed, report.
PyObject *selftest(PyObject *const *) {
  who = "selftest";
  static char sink[1 << 16];
  int nf = 0;
  const int n = clab_selftest(sink, static_cast<int>(sizeof sink), &nf);
  Obj out = dict();
  put(out, "n", own(PyLong_FromLong(n)));
  put(out, "nfail", own(PyLong_FromLong(nf)));
  put(out, "report", own(PyUnicode_FromString(sink)));
  return out.release();
}

// as.character() of one number, for labels made in Python.
PyObject *rchar(PyObject *const *a) {
  who = "rchar";
  const double x = PyFloat_AsDouble(a[0]);
  if (x == -1.0 && PyErr_Occurred()) throw Raised{};
  return PyUnicode_FromString(rfmt(x).c_str());
}

// A hook for the tests, on lengths no test can allocate: with a name, the
// length check of the binding; with False, that of the kernels; with None,
// the length as a message prints it.
PyObject *tlen(PyObject *const *a) {
  who = "clab_tlen";
  const double v = PyFloat_AsDouble(a[0]);
  if (v == -1.0 && PyErr_Occurred()) throw Raised{};
  if (a[1] == Py_None) return PyUnicode_FromString(rfmt(v).c_str());
  const char *s = PyUnicode_Check(a[1]) ? PyUnicode_AsUTF8AndSize(a[1], nullptr) : nullptr;
  if (s == nullptr && PyErr_Occurred()) throw Raised{};
  return PyLong_FromLong(ilen(v, s));
}

using Fn = PyObject *(*)(PyObject *const *);

// Entry point f, which takes n arguments, called with the tuple args: the
// one place an exception becomes a Python error.
PyObject *call(Fn f, int n, PyObject *args) {
  if (PyTuple_Size(args) != n) {
    PyErr_Format(PyExc_TypeError, "expected %d arguments, got %zd", n, PyTuple_Size(args));
    return nullptr;
  }
  PyObject *a[8] = {nullptr};
  for (int i = 0; i < n; ++i) a[i] = PyTuple_GetItem(args, i);
  try {
    return f(a);
  } catch (const Raised &) {
    return nullptr;
  } catch (const std::exception &) {  // std::bad_alloc, std::length_error
    return PyErr_NoMemory();
  }
}

// Adding a function: its entry point above, one line here, and its Python
// wrapper in crimlab/.
#define CLAB_FN(f, n) \
  {#f, [](PyObject *, PyObject *a) { return call(f, n, a); }, METH_VARARGS, nullptr}
PyMethodDef methods[] = {
    CLAB_FN(vodtst, 5),
    CLAB_FN(orrisk, 3),
    CLAB_FN(outtst, 3),
    CLAB_FN(dspbm, 5),
    CLAB_FN(dyadrt, 5),
    CLAB_FN(gtmp, 5),
    CLAB_FN(cllr, 2),
    CLAB_FN(lrpost, 2),
    CLAB_FN(mixlr, 6),
    CLAB_FN(dbsrch, 3),
    CLAB_FN(cllrmin, 2),
    CLAB_FN(lrcal, 2),
    CLAB_FN(tipplt, 2),
    CLAB_FN(eceplt, 3),
    CLAB_FN(kinidx, 5),
    CLAB_FN(patidx, 5),
    CLAB_FN(clrsmp, 5),
    CLAB_FN(expbnd, 5),
    CLAB_FN(misrd, 5),
    CLAB_FN(rtmhs, 4),
    CLAB_FN(incap, 5),
    CLAB_FN(incapr, 4),
    CLAB_FN(sellab, 4),
    CLAB_FN(losbnd, 2),
    CLAB_FN(sntbnd, 5),
    CLAB_FN(hzdsel, 5),
    CLAB_FN(selftest, 0),
    CLAB_FN(rchar, 1),
    CLAB_FN(tlen, 2),
    {nullptr, nullptr, 0, nullptr}};
#undef CLAB_FN

PyModuleDef module = {PyModuleDef_HEAD_INIT, "_clab", "The compiled core of crimlab.", 0,
                      methods, nullptr, nullptr, nullptr, nullptr};

}  // namespace

PyMODINIT_FUNC PyInit__clab(void) { return PyModuleDef_Init(&module); }
