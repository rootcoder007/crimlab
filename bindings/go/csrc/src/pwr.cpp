// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The readers of hwpwr() that work from files and programs (src/pwr.h):
// the Linux powercap (RAPL), DRM and power-supply trees, the kernel's CPU
// accounting, the NVIDIA management library loaded at run time with
// nvidia-smi as the last resort, and the snapshot text. The macOS and
// Windows readers are in pwrmac.cpp and pwrwin.cpp.
//
// Units, as the kernel documents them (Documentation/power/powercap,
// Documentation/ABI/testing/sysfs-class-power, the amdgpu and Intel GPU
// hwmon pages): energy_uj and energy1_input in microjoules,
// max_energy_range_uj the counter's range, power*_ files in microwatts,
// current_now in microamperes and voltage_now in microvolts.
#ifdef _WIN32
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#include "pwr.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
// WebAssembly (WASI) has no processes and no shared libraries: there no
// program is run and no library is loaded, and each reader that needs one
// reports that it is not available
#ifndef __wasi__
#include <dlfcn.h>
#include <sys/wait.h>
#endif
#endif

namespace pwr {

bool Cfg::on(const char *b) const { return std::find(skip.begin(), skip.end(), b) == skip.end(); }

// ---- files, folders, programs --------------------------------------------

bool rd(const std::string &path, std::string &out, std::string &err) {
  out.clear();
  err.clear();
  errno = 0;
  std::FILE *f = std::fopen(path.c_str(), "rb");
  if (f == nullptr) {
    err = errno == EACCES || errno == EPERM ? "denied" : errno == ENOENT ? "absent" : std::strerror(errno);
    return false;
  }
  char buf[4096];
  std::size_t k;
  while ((k = std::fread(buf, 1, sizeof buf, f)) > 0) out.append(buf, k);
  // a sysfs attribute can refuse at read() rather than at open()
  const bool ok = !std::ferror(f);
  if (!ok) err = errno == EACCES || errno == EPERM ? "denied" : std::strerror(errno);
  std::fclose(f);
  return ok;
}

std::vector<std::string> ls(const std::string &dir) {
  std::vector<std::string> v;
#ifdef _WIN32
  WIN32_FIND_DATAA f;
  HANDLE h = FindFirstFileA((dir + "\\*").c_str(), &f);
  if (h == INVALID_HANDLE_VALUE) return v;
  do {
    if (std::strcmp(f.cFileName, ".") != 0 && std::strcmp(f.cFileName, "..") != 0) v.emplace_back(f.cFileName);
  } while (FindNextFileA(h, &f));
  FindClose(h);
#else
  DIR *d = opendir(dir.c_str());
  if (d == nullptr) return v;
  while (dirent *e = readdir(d)) {
    if (std::strcmp(e->d_name, ".") != 0 && std::strcmp(e->d_name, "..") != 0) v.emplace_back(e->d_name);
  }
  closedir(d);
#endif
  std::sort(v.begin(), v.end());
  return v;
}

bool run(const std::string &cmd, std::string &out, int &st) {
  out.clear();
  st = -1;
#ifdef __wasi__
  (void)cmd;
  return false;
#else
#ifdef _WIN32
  std::FILE *p = _popen((cmd + " 2>NUL").c_str(), "r");
#else
  std::FILE *p = popen((cmd + " 2>/dev/null").c_str(), "r");
#endif
  if (p == nullptr) return false;
  char buf[4096];
  std::size_t k;
  while ((k = std::fread(buf, 1, sizeof buf, p)) > 0) out.append(buf, k);
#ifdef _WIN32
  st = _pclose(p);
#else
  const int w = pclose(p);
  st = w == -1 ? -1 : WIFEXITED(w) ? WEXITSTATUS(w) : -1;
  // the shell's "command not found"
  if (st == 127) st = -1;
#endif
  return st == 0;
#endif
}

namespace {

std::string trim(const std::string &s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

// A whole decimal number (the kernel's format for counters), or false.
bool num(const std::string &s, double &v) {
  const std::string t = trim(s);
  if (t.empty()) return false;
  char *e = nullptr;
  errno = 0;
  v = std::strtod(t.c_str(), &e);
  return e != nullptr && *e == '\0' && errno == 0 && std::isfinite(v);
}

// The first line of a file, trimmed; "" when it cannot be read.
std::string line1(const std::string &p) {
  std::string s, e;
  if (!rd(p, s, e)) return "";
  return trim(s.substr(0, s.find('\n')));
}

bool fnum(const std::string &p, double &v, std::string &err) {
  std::string s;
  if (!rd(p, s, err)) return false;
  if (!num(s, v)) {
    err = "malformed";
    return false;
  }
  return true;
}

bool isdig(const std::string &s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) { return c >= '0' && c <= '9'; });
}

// A program path quoted for the shell.
std::string q(const std::string &p) {
#ifdef _WIN32
  return "\"" + p + "\"";
#else
  std::string r = "'";
  for (char c : p) r += c == '\'' ? std::string("'\\''") : std::string(1, c);
  return r + "'";
#endif
}

std::string clean(std::string s) {
  for (char &c : s) {
    if (c == '\t' || c == '\n' || c == '\r') c = ' ';
  }
  return s;
}

}  // namespace

// ---- Linux: RAPL through powercap ----------------------------------------

void rapl(const Cfg &c, Snap &s) {
  const std::string base = c.root + "/sys/class/powercap";
  // the CPU vendor, from the kernel's description of the processor
  std::string ci, e, ven;
  if (rd(c.root + "/proc/cpuinfo", ci, e)) {
    if (ci.find("GenuineIntel") != std::string::npos) ven = "intel";
    if (ci.find("AuthenticAMD") != std::string::npos || ci.find("HygonGenuine") != std::string::npos) ven = "amd";
  }
  // Top-level zones "intel-rapl:N" (and the MMIO copy "intel-rapl-mmio:N")
  // whose name starts with "package" are the CPU packages; their subzone
  // named "dram" is memory. core, uncore and psys are never added: core
  // and uncore are parts of a package, psys covers more than the CPU.
  struct Cand {
    std::string path, key, name;
    bool mmio, mem;
  };
  std::vector<Cand> cand;
  bool any = false;
  for (const std::string &z : ls(base)) {
    const bool mmio = z.rfind("intel-rapl-mmio:", 0) == 0;
    if (!mmio && z.rfind("intel-rapl:", 0) != 0) continue;
    any = true;
    const std::string rest = z.substr(z.find(':') + 1);
    const std::size_t k = rest.find(':');
    const std::string nm = line1(base + "/" + z + "/name");
    if (k == std::string::npos && isdig(rest)) {
      if (nm.rfind("package", 0) == 0) cand.push_back({base + "/" + z, nm, nm, mmio, false});
    } else if (k != std::string::npos && isdig(rest.substr(0, k)) && isdig(rest.substr(k + 1)) && nm == "dram") {
      const std::string par = line1(base + "/" + z.substr(0, z.rfind(':')) + "/name");
      cand.push_back({base + "/" + z, "dram:" + par, par.empty() ? "dram" : "dram of " + par, mmio, true});
    }
  }
  if (cand.empty()) {
    Dev d;
    d.id = "rapl:none";
    d.kind = "cpu";
    d.vendor = ven;
    d.name = "CPU";
    d.bad = 1;
    d.why = any ? "not measurable here: the powercap tree has no CPU package domain"
                : "not measurable here: no RAPL powercap domains (the case in most virtual machines)";
    s.d.push_back(d);
    return;
  }
  // One zone per key: the same package can appear through MSR and MMIO;
  // a readable copy first, then the MMIO one.
  std::map<std::string, Dev> got;
  std::map<std::string, std::pair<bool, bool>> rank;  // key -> (readable, mmio)
  bool dram = false;
  for (const Cand &z : cand) {
    Dev d;
    d.id = "rapl:" + z.key;
    d.kind = z.mem ? "memory" : "cpu";
    d.vendor = ven;
    d.name = z.name;
    d.aon = 0;
    double v = 0, mx = 0;
    std::string err;
    if (fnum(z.path + "/energy_uj", v, err) && v >= 0) {
      d.src = 1;
      d.ctr = v * 1e-6;
      std::string e2;
      if (fnum(z.path + "/max_energy_range_uj", mx, e2) && mx > 0) d.wrap = mx * 1e-6;
    } else {
      d.bad = 1;
      d.why = err == "denied"   ? "not readable: needs read access to the powercap files"
              : err == "absent" ? "not measurable here: the zone has no energy counter"
                                : "unreadable energy counter (" + (err.empty() ? std::string("malformed") : err) + ")";
    }
    if (!z.mem) {
      // the package's rated power as the machine reports it: the maximum of
      // the long-term constraint, else its limit
      double w = 0;
      std::string e3;
      if (fnum(z.path + "/constraint_0_max_power_uw", w, e3) && w > 0) {
        d.rated = w * 1e-6;
      } else if (line1(z.path + "/constraint_0_name") == "long_term" &&
                 fnum(z.path + "/constraint_0_power_limit_uw", w, e3) && w > 0) {
        d.rated = w * 1e-6;
      }
    } else {
      dram = true;
    }
    const std::pair<bool, bool> r(d.src == 1, z.mmio);
    auto it = rank.find(z.key);
    if (it == rank.end() || r > it->second) {
      rank[z.key] = r;
      got[z.key] = d;
    }
  }
  for (auto &kv : got) s.d.push_back(kv.second);
  if (!dram) {
    Dev d;
    d.id = "rapl:dram";
    d.kind = "memory";
    d.vendor = ven;
    d.name = "memory";
    d.bad = 1;
    d.why = "not measurable here: the CPU exposes no DRAM energy domain";
    s.d.push_back(d);
  }
}

// ---- Linux: GPUs through DRM ---------------------------------------------

void drm(const Cfg &c, Snap &s) {
  const std::string base = c.root + "/sys/class/drm";
  for (const std::string &cd : ls(base)) {
    if (cd.rfind("card", 0) != 0 || !isdig(cd.substr(4))) continue;
    const std::string dev = base + "/" + cd + "/device";
    const std::string vid = line1(dev + "/vendor"), pid = line1(dev + "/device");
    if (vid.empty()) continue;
    Dev d;
    d.id = "drm:" + cd;
    d.kind = "gpu";
    std::string pci = vid.size() > 2 ? vid.substr(2) : vid;
    if (pid.size() > 2) pci += ":" + pid.substr(2);
    d.name = cd + " (" + pci + ")";
    // the hwmon folder of the device, and its name
    std::string hw, hwn;
    const std::vector<std::string> hws = ls(dev + "/hwmon");
    if (!hws.empty()) {
      hw = dev + "/hwmon/" + hws[0];
      hwn = line1(hw + "/name");
    }
    double v = 0;
    std::string err;
    if (vid == "0x1002") {
      d.vendor = "amd";
      if (fnum(dev + "/gpu_busy_percent", v, err) && v >= 0 && v <= 100) d.util = v / 100;
      if (!hw.empty() && hwn == "amdgpu" &&
          ((fnum(hw + "/power1_average", v, err) && v >= 0) || (fnum(hw + "/power1_input", v, err) && v >= 0))) {
        d.src = 2;
        d.w = v * 1e-6;
        d.why = "on an APU this figure includes the CPU";
      } else {
        d.bad = 1;
        d.why = "present, power unknown: the driver exposes no power file";
      }
    } else if (vid == "0x8086") {
      d.vendor = "intel";
      if (!hw.empty() && fnum(hw + "/energy1_input", v, err) && v >= 0) {
        d.src = 1;
        d.ctr = v * 1e-6;
      } else {
        d.bad = 1;
        // Intel's integrated graphics sit at PCI device 2, function 0, of bus 0
        std::string ue, e2;
        rd(dev + "/uevent", ue, e2);
        const std::size_t k = ue.find("PCI_SLOT_NAME=");
        const std::string slot = k == std::string::npos ? "" : trim(ue.substr(k + 14, ue.find('\n', k) - k - 14));
        if (slot.size() >= 7 && slot.compare(slot.size() - 7, 7, "00:02.0") == 0) {
          d.inpkg = 1;
          d.why = "integrated: its energy is part of the CPU package figure";
        } else {
          d.why = "present, power unknown: the kernel exposes no energy file for it";
        }
      }
    } else if (vid == "0x10de") {
      d.vendor = "nvidia";
      d.bad = 1;
      d.id = "drmnv:" + cd;
      d.why = "present, power unknown: the NVIDIA driver library was not found or gave no reading";
    } else {
      d.bad = 1;
      d.why = "present, power unknown: no reader for this vendor";
    }
    s.d.push_back(d);
  }
}

// ---- Linux: batteries ------------------------------------------------------

void bat(const Cfg &c, Snap &s) {
  const std::string base = c.root + "/sys/class/power_supply";
  for (const std::string &b : ls(base)) {
    const std::string p = base + "/" + b;
    if (line1(p + "/type") != "Battery") continue;
    Dev d;
    d.id = "bat:" + b;
    d.kind = "machine";
    d.name = "battery " + b;
    const std::string st = line1(p + "/status");
    double pw = 0, cu = 0, vo = 0;
    std::string e;
    if (st != "Discharging") {
      d.bad = 1;
      d.why = "not on battery (status " + (st.empty() ? std::string("unknown") : st) +
              "): the battery gives the machine's power only while it discharges";
    } else if (fnum(p + "/power_now", pw, e) && pw > 0) {
      d.src = 2;
      d.w = pw * 1e-6;
    } else if (fnum(p + "/current_now", cu, e) && fnum(p + "/voltage_now", vo, e) && cu > 0 && vo > 0) {
      d.src = 2;
      d.w = cu * vo * 1e-12;
    } else {
      d.bad = 1;
      d.why = "the battery reports no discharge power";
    }
    s.d.push_back(d);
  }
}

// ---- Linux: CPU time of the system and of a process tree -----------------

namespace {

// Fields of /proc/<pid>/stat after the command name: state is field 3, so
// field k (1-based, as proc(5) numbers them) is f[k - 3].
bool pstat(const std::string &p, std::vector<std::string> &f) {
  std::string s, e;
  if (!rd(p, s, e)) return false;
  const std::size_t r = s.rfind(')');
  if (r == std::string::npos) return false;
  f.clear();
  std::size_t i = r + 1;
  while (i < s.size()) {
    while (i < s.size() && s[i] == ' ') ++i;
    std::size_t j = i;
    while (j < s.size() && s[j] != ' ' && s[j] != '\n') ++j;
    if (j > i) f.push_back(s.substr(i, j - i));
    if (j < s.size() && s[j] == '\n') break;
    i = j;
  }
  return f.size() >= 15;
}

}  // namespace

void proc(const Cfg &c, Snap &s) {
  std::string st, e;
  if (rd(c.root + "/proc/stat", st, e) && st.rfind("cpu ", 0) == 0) {
    // user nice system idle iowait irq softirq steal guest guest_nice; guest
    // time is already inside user and nice
    std::vector<double> v;
    const std::string l = st.substr(4, st.find('\n') - 4);
    char *p = const_cast<char *>(l.c_str());
    for (int k = 0; k < 8; ++k) {
      char *q2 = nullptr;
      const double x = std::strtod(p, &q2);
      if (q2 == p) break;
      v.push_back(x);
      p = q2;
    }
    if (v.size() >= 4) {
      double busy = v[0] + v[1] + v[2];
      for (std::size_t k = 5; k < v.size(); ++k) busy += v[k];
      s.scpu = busy;
      s.stot = busy + v[3] + (v.size() > 4 ? v[4] : 0);
    }
  }
  if (c.pid <= 0) return;
  // the process and every live descendant; reaped children are in cutime
  // and cstime of their parent
  std::map<long, long> par;
  std::map<long, double> tm;
  std::vector<std::string> f;
  for (const std::string &d : ls(c.root + "/proc")) {
    if (!isdig(d) || !pstat(c.root + "/proc/" + d + "/stat", f)) continue;
    const long id = std::atol(d.c_str());
    par[id] = std::atol(f[1].c_str());
    tm[id] = std::atof(f[11].c_str()) + std::atof(f[12].c_str()) + std::atof(f[13].c_str()) + std::atof(f[14].c_str());
  }
  if (!tm.count(c.pid)) return;
  std::set<long> in{c.pid};
  for (bool more = true; more;) {
    more = false;
    for (auto &kv : par) {
      if (!in.count(kv.first) && in.count(kv.second)) {
        in.insert(kv.first);
        more = true;
      }
    }
  }
  double t = 0;
  for (long id : in) t += tm[id];
  s.pcpu = t;
}

// ---- NVIDIA ----------------------------------------------------------------

namespace {

void *lopen(const char *p) {
#if defined(_WIN32)
  return reinterpret_cast<void *>(LoadLibraryA(p));
#elif defined(__wasi__)
  (void)p;
  return nullptr;
#else
  return dlopen(p, RTLD_NOW | RTLD_LOCAL);
#endif
}

void *lsym(void *h, const char *n) {
#if defined(_WIN32)
  return reinterpret_cast<void *>(GetProcAddress(static_cast<HMODULE>(h), n));
#elif defined(__wasi__)
  (void)h;
  (void)n;
  return nullptr;
#else
  return dlsym(h, n);
#endif
}

// NVML's types as its reference manual gives them; nvmlDevice_t is opaque.
typedef void *nvdev;
struct nvutl {
  unsigned int gpu, memory;
};
struct nvpus {
  unsigned int pid;
  unsigned long long timeStamp;
  unsigned int smUtil, memUtil, encUtil, decUtil;
};
typedef int (*f_v)();
typedef int (*f_u)(unsigned int *);
typedef int (*f_h)(unsigned int, nvdev *);
typedef int (*f_n)(nvdev, char *, unsigned int);
typedef int (*f_e)(nvdev, unsigned long long *);
typedef int (*f_p)(nvdev, unsigned int *);
typedef int (*f_r)(nvdev, nvutl *);
typedef int (*f_s)(nvdev, nvpus *, unsigned int *, unsigned long long);

struct Nvml {
  bool tried = false, up = false;
  std::string lib, why;
  f_u count = nullptr;
  f_h handle = nullptr;
  f_n name = nullptr;
  f_e energy = nullptr;
  f_p power = nullptr;
  f_r util = nullptr;
  f_s putil = nullptr;
};

// The library is loaded and started once per process and never unloaded:
// the driver keeps its state for the next reading.
Nvml &nvml(const std::string &want) {
  static std::map<std::string, Nvml> st;
  Nvml &n = st[want];
  if (n.tried) return n;
  n.tried = true;
  std::vector<std::string> names;
  if (!want.empty()) {
    names.push_back(want);
  } else {
#ifdef _WIN32
    names = {"nvml.dll"};
    const char *pf = std::getenv("ProgramFiles");
    if (pf != nullptr) names.push_back(std::string(pf) + "\\NVIDIA Corporation\\NVSMI\\nvml.dll");
#else
    names = {"libnvidia-ml.so.1", "libnvidia-ml.so"};
#endif
  }
  void *h = nullptr;
  for (const std::string &p : names) {
    if ((h = lopen(p.c_str())) != nullptr) {
      n.lib = p;
      break;
    }
  }
  if (h == nullptr) return n;
  // every symbol is optional; without init, count and handle nothing is read
  f_v init = reinterpret_cast<f_v>(lsym(h, "nvmlInit_v2"));
  if (init == nullptr) init = reinterpret_cast<f_v>(lsym(h, "nvmlInit"));
  n.count = reinterpret_cast<f_u>(lsym(h, "nvmlDeviceGetCount_v2"));
  if (n.count == nullptr) n.count = reinterpret_cast<f_u>(lsym(h, "nvmlDeviceGetCount"));
  n.handle = reinterpret_cast<f_h>(lsym(h, "nvmlDeviceGetHandleByIndex_v2"));
  if (n.handle == nullptr) n.handle = reinterpret_cast<f_h>(lsym(h, "nvmlDeviceGetHandleByIndex"));
  n.name = reinterpret_cast<f_n>(lsym(h, "nvmlDeviceGetName"));
  n.energy = reinterpret_cast<f_e>(lsym(h, "nvmlDeviceGetTotalEnergyConsumption"));
  n.power = reinterpret_cast<f_p>(lsym(h, "nvmlDeviceGetPowerUsage"));
  n.util = reinterpret_cast<f_r>(lsym(h, "nvmlDeviceGetUtilizationRates"));
  n.putil = reinterpret_cast<f_s>(lsym(h, "nvmlDeviceGetProcessUtilization"));
  if (init == nullptr || n.count == nullptr || n.handle == nullptr) {
    n.why = "the NVIDIA driver library lacks the functions to enumerate GPUs";
    return n;
  }
  const int rc = init();
  if (rc != 0) {
    n.why = "the NVIDIA driver library did not start (code " + std::to_string(rc) + ")";
    return n;
  }
  n.up = true;
  return n;
}

}  // namespace

void smi(const std::string &out, Snap &s) {
  std::size_t i = 0;
  int k = 0;
  while (i < out.size()) {
    std::size_t j = out.find('\n', i);
    if (j == std::string::npos) j = out.size();
    const std::string l = trim(out.substr(i, j - i));
    i = j + 1;
    if (l.empty()) continue;
    const std::size_t a = l.find(','), b = a == std::string::npos ? a : l.find(',', a + 1);
    Dev d;
    d.id = "smi:" + std::to_string(k++);
    d.kind = "gpu";
    d.vendor = "nvidia";
    if (b == std::string::npos) {
      d.name = "NVIDIA GPU";
      d.bad = 1;
      d.why = "nvidia-smi gave a line that is not power, utilisation, name";
      s.d.push_back(d);
      continue;
    }
    d.name = clean(trim(l.substr(b + 1)));
    double w = 0, u = 0;
    if (num(l.substr(0, a), w) && w >= 0) {
      d.src = 2;
      d.w = w;
    } else {
      d.bad = 1;
      d.why = "nvidia-smi reports no power for this GPU";
    }
    if (num(l.substr(a + 1, b - a - 1), u) && u >= 0 && u <= 100) d.util = u / 100;
    s.d.push_back(d);
  }
}

void nv(const Cfg &c, Snap &s) {
  static Nvml none;
  Nvml &n = c.on("nvml") ? nvml(c.nvml) : none;
  if (n.up) {
    unsigned int cnt = 0;
    if (n.count(&cnt) == 0) {
      const auto now = std::chrono::system_clock::now().time_since_epoch();
      const unsigned long long us =
          static_cast<unsigned long long>(std::chrono::duration_cast<std::chrono::microseconds>(now).count());
      for (unsigned int i = 0; i < cnt; ++i) {
        nvdev h = nullptr;
        if (n.handle(i, &h) != 0 || h == nullptr) continue;
        Dev d;
        d.id = "nvml:" + std::to_string(i);
        d.kind = "gpu";
        d.vendor = "nvidia";
        char nm[256] = {0};
        d.name = n.name != nullptr && n.name(h, nm, sizeof nm - 1) == 0 ? clean(nm) : "NVIDIA GPU " + std::to_string(i);
        unsigned long long mj = 0;
        unsigned int mw = 0;
        if (n.energy != nullptr && n.energy(h, &mj) == 0) {
          d.src = 1;
          d.ctr = static_cast<double>(mj) * 1e-3;
        }
        if (n.power != nullptr && n.power(h, &mw) == 0) {
          d.w = mw * 1e-3;
          if (d.src == 0) d.src = 2;
        }
        if (d.src == 0) {
          d.bad = 1;
          d.why = "the driver reports no energy or power for this GPU";
        }
        nvutl u{};
        if (n.util != nullptr && n.util(h, &u) == 0 && u.gpu <= 100) d.util = u.gpu / 100.0;
        // the process's share of the SM samples of the last 10 seconds
        if (n.putil != nullptr && c.pid > 0) {
          unsigned int m = 0;
          const unsigned long long since = us > 10000000ULL ? us - 10000000ULL : 0;
          n.putil(h, nullptr, &m, since);
          if (m > 0 && m < 100000) {
            std::vector<nvpus> v(m);
            if (n.putil(h, v.data(), &m, since) == 0) {
              double mine = 0, all = 0;
              for (unsigned int k = 0; k < m && k < v.size(); ++k) {
                all += v[k].smUtil;
                if (static_cast<long>(v[k].pid) == c.pid) mine += v[k].smUtil;
              }
              if (all > 0) d.pown = mine / all;
            }
          }
        }
        s.d.push_back(d);
      }
      return;
    }
  }
  if (!n.why.empty()) {
    Dev d;
    d.id = "nvml:lib";
    d.kind = "gpu";
    d.vendor = "nvidia";
    d.name = "NVIDIA GPU";
    d.bad = 1;
    d.why = n.why;
    s.d.push_back(d);
    return;
  }
  // the last resort: nvidia-smi, when the library was not found
  if (!c.on("smi") || c.smi.empty()) return;
  std::string out;
  int st = 0;
  const bool ok = run(q(c.smi) + " --query-gpu=power.draw,utilization.gpu,name --format=csv,noheader,nounits", out, st);
  if (ok) {
    smi(out, s);
  } else if (st > 0) {
    Dev d;
    d.id = "smi:fail";
    d.kind = "gpu";
    d.vendor = "nvidia";
    d.name = "NVIDIA GPU";
    d.bad = 1;
    d.why = "nvidia-smi failed (exit status " + std::to_string(st) + ")";
    s.d.push_back(d);
  }
}

// ---- parsers of ioreg (macOS) and of the EMI metadata (Windows) ----------

std::vector<double> ioacc(const std::string &out) {
  std::vector<double> v;
  const std::string key = "\"Device Utilization %\"=";
  for (std::size_t i = out.find(key); i != std::string::npos; i = out.find(key, i + 1)) {
    const char *p = out.c_str() + i + key.size();
    char *e = nullptr;
    const double x = std::strtod(p, &e);
    if (e != p && x >= 0 && x <= 100) v.push_back(x / 100);
  }
  return v;
}

double iobat(const std::string &out, std::string &why) {
  // top-level keys of the battery's registry entry: one per line, "key" = value
  std::map<std::string, std::string> kv;
  std::size_t i = 0;
  while (i < out.size()) {
    std::size_t j = out.find('\n', i);
    if (j == std::string::npos) j = out.size();
    const std::string l = trim(out.substr(i, j - i));
    i = j + 1;
    const std::size_t e = l.find("\" = ");
    if (l.size() > 1 && l[0] == '"' && e != std::string::npos) kv[l.substr(1, e - 1)] = trim(l.substr(e + 4));
  }
  if (!kv.count("Voltage") || (!kv.count("InstantAmperage") && !kv.count("Amperage"))) {
    why = "no battery found";
    return -1;
  }
  if (kv["ExternalConnected"] == "Yes") {
    why = "on mains power: the battery gives the machine's power only while it discharges";
    return -1;
  }
  // a negative current is printed as its 64-bit two's complement
  const std::string as = kv.count("InstantAmperage") ? kv["InstantAmperage"] : kv["Amperage"];
  double v = 0;
  const std::string at = trim(as);
  char *e = nullptr;
  errno = 0;
  const unsigned long long ua = at.empty() || at[0] == '-' ? 0ULL : std::strtoull(at.c_str(), &e, 10);
  if (at.empty() || at[0] == '-' || e == nullptr || *e != '\0' || errno != 0 || !num(kv["Voltage"], v) || v <= 0) {
    why = "the battery's current or voltage is malformed";
    return -1;
  }
  const double a = static_cast<double>(static_cast<std::int64_t>(ua));
  if (a > 0) {
    why = "the battery is charging";
    return -1;
  }
  why.clear();
  return -a * v * 1e-6;
}

namespace {

std::string u16(const unsigned char *b, unsigned long nbytes) {
  std::string s;
  for (unsigned long i = 0; i + 1 < nbytes; i += 2) {
    const unsigned c = b[i] | (b[i + 1] << 8);
    if (c == 0) break;
    s += c < 128 ? static_cast<char>(c) : '?';
  }
  return s;
}

unsigned rd16(const unsigned char *b) { return b[0] | (b[1] << 8); }
unsigned long rd32(const unsigned char *b) {
  return b[0] | (b[1] << 8) | (static_cast<unsigned long>(b[2]) << 16) | (static_cast<unsigned long>(b[3]) << 24);
}

}  // namespace

bool emimeta(const unsigned char *b, unsigned long n, int ver, std::vector<std::string> &names,
             std::vector<int> &pwh) {
  names.clear();
  pwh.clear();
  // EMI_METADATA_V1: unit (4), HardwareOEM and HardwareModel (16 WCHAR
  // each), HardwareRevision, MeteredHardwareNameSize, then the name.
  if (ver == 1) {
    if (n < 72) return false;
    const unsigned long sz = rd16(b + 70);
    if (72 + sz > n) return false;
    names.push_back(u16(b + 72, sz));
    pwh.push_back(rd32(b) == 0);
    return true;
  }
  // EMI_METADATA_V2: HardwareOEM, HardwareModel, HardwareRevision,
  // ChannelCount, then the channels, each unit (4), ChannelNameSize (2)
  // and the name, packed one after the other.
  if (ver != 2 || n < 68) return false;
  const unsigned cnt = rd16(b + 66);
  unsigned long at = 68;
  for (unsigned k = 0; k < cnt; ++k) {
    if (at + 6 > n) return false;
    const unsigned long sz = rd16(b + at + 4);
    if (at + 6 + sz > n) return false;
    names.push_back(u16(b + at + 6, sz));
    pwh.push_back(rd32(b + at) == 0);
    at += 6 + sz;
  }
  return true;
}

// ---- the snapshot ----------------------------------------------------------

Snap snap(const Cfg &c) {
  Snap s;
  s.t = std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
#ifdef __linux__
  const bool lx = true;
#else
  const bool lx = !c.root.empty();
#endif
  if (lx) {
    if (c.on("rapl")) rapl(c, s);
    if (c.on("gpu")) drm(c, s);
    if (c.on("battery")) bat(c, s);
    if (c.on("proc")) proc(c, s);
  }
  nv(c, s);
  mac(c, s);
  win(c, s);
  // a DRM card of NVIDIA's stays only when no NVIDIA reader answered
  bool nvd = false;
  for (const Dev &d : s.d) nvd = nvd || (d.vendor == "nvidia" && d.id.rfind("drmnv:", 0) != 0 && !d.bad);
  if (nvd) {
    s.d.erase(std::remove_if(s.d.begin(), s.d.end(), [](const Dev &d) { return d.id.rfind("drmnv:", 0) == 0; }),
              s.d.end());
  }
  return s;
}

std::string put(const Snap &s) {
  char b[64];
  auto g = [&](double x) {
    std::snprintf(b, sizeof b, "%.17g", x);
    return std::string(b);
  };
  std::string o = "hwpwr1\t" + g(s.t) + "\t" + g(s.pcpu) + "\t" + g(s.scpu) + "\t" + g(s.stot) + "\n";
  for (const Dev &d : s.d) {
    o += clean(d.id) + "\t" + clean(d.kind) + "\t" + clean(d.vendor) + "\t" + clean(d.name) + "\t" + clean(d.why) +
         "\t" + std::to_string(d.src) + "\t" + g(d.ctr) + "\t" + g(d.wrap) + "\t" + g(d.w) + "\t" + g(d.util) + "\t" +
         g(d.rated) + "\t" + g(d.pown) + "\t" + std::to_string(d.inpkg) + "\t" + std::to_string(d.aon) + "\t" +
         std::to_string(d.bad) + "\n";
  }
  return o;
}

bool get(const std::string &txt, Snap &s) {
  s = Snap();
  std::vector<std::vector<std::string>> ln;
  std::size_t i = 0;
  while (i < txt.size()) {
    std::size_t j = txt.find('\n', i);
    if (j == std::string::npos) j = txt.size();
    std::vector<std::string> f;
    std::size_t a = i;
    for (std::size_t k = i; k <= j; ++k) {
      if (k == j || txt[k] == '\t') {
        f.push_back(txt.substr(a, k - a));
        a = k + 1;
      }
    }
    if (j > i) ln.push_back(f);
    i = j + 1;
  }
  if (ln.empty() || ln[0].size() != 5 || ln[0][0] != "hwpwr1") return false;
  double h[4];
  for (int k = 0; k < 4; ++k) {
    if (!num(ln[0][k + 1], h[k])) return false;
  }
  s.t = h[0];
  s.pcpu = h[1];
  s.scpu = h[2];
  s.stot = h[3];
  for (std::size_t r = 1; r < ln.size(); ++r) {
    const std::vector<std::string> &f = ln[r];
    if (f.size() != 15) return false;
    Dev d;
    d.id = f[0];
    d.kind = f[1];
    d.vendor = f[2];
    d.name = f[3];
    d.why = f[4];
    double x[10];
    for (int k = 0; k < 10; ++k) {
      if (!num(f[k + 5], x[k])) return false;
    }
    d.src = static_cast<int>(x[0]);
    d.ctr = x[1];
    d.wrap = x[2];
    d.w = x[3];
    d.util = x[4];
    d.rated = x[5];
    d.pown = x[6];
    d.inpkg = static_cast<int>(x[7]);
    d.aon = static_cast<int>(x[8]);
    d.bad = static_cast<int>(x[9]);
    if (d.id.empty() || d.src < 0 || d.src > 2) return false;
    s.d.push_back(d);
  }
  return true;
}

}  // namespace pwr
