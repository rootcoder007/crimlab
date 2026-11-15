// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The macOS readers of hwpwr() (src/pwr.h). On Apple silicon the energy
// counters of the IOReport group "Energy Model" (CPU, GPU, neural engine,
// DRAM) are read without administrator rights from libIOReport, loaded at
// run time with CoreFoundation, so nothing is linked. IOReport is not a
// documented public interface of Apple's: every symbol is optional, every
// object is checked and released, and any failure leaves the counters out
// with the reason. The GPU's busy share comes from ioreg (the
// IOAccelerator's "Device Utilization %"), the battery's discharge from
// ioreg's AppleSmartBattery entry (current times voltage), and the CPU time
// of the process and of the system from getrusage() and the host's tick
// counts.
#include "pwr.h"

#ifdef __APPLE__
#include <CoreFoundation/CoreFoundation.h>
#include <dlfcn.h>
#include <mach/mach.h>
#include <sys/resource.h>
#include <sys/sysctl.h>
#include <unistd.h>

#include <cstring>
#include <map>
#endif

namespace pwr {

#ifndef __APPLE__

void mac(const Cfg &, Snap &) {}

#else

namespace {

typedef CFStringRef (*cf_mks)(CFAllocatorRef, const char *, CFStringEncoding);
typedef Boolean (*cf_gets)(CFStringRef, char *, CFIndex, CFStringEncoding);
typedef void (*cf_rel)(CFTypeRef);
typedef const void *(*cf_dget)(CFDictionaryRef, const void *);
typedef CFIndex (*cf_dcnt)(CFDictionaryRef);
typedef CFMutableDictionaryRef (*cf_dcpy)(CFAllocatorRef, CFIndex, CFDictionaryRef);
typedef CFIndex (*cf_acnt)(CFArrayRef);
typedef const void *(*cf_aget)(CFArrayRef, CFIndex);
typedef CFTypeID (*cf_tid)(CFTypeRef);
typedef CFTypeID (*cf_atid)();

typedef CFDictionaryRef (*io_grp)(CFStringRef, CFStringRef, uint64_t, uint64_t, uint64_t);
typedef void *(*io_sub)(void *, CFMutableDictionaryRef, CFMutableDictionaryRef *, uint64_t, CFTypeRef);
typedef CFDictionaryRef (*io_smp)(void *, CFMutableDictionaryRef, CFTypeRef);
typedef int64_t (*io_int)(CFDictionaryRef, int32_t);
typedef CFStringRef (*io_str)(CFDictionaryRef);

struct Ior {
  bool tried = false, up = false;
  std::string why;
  cf_mks mks = nullptr;
  cf_gets gets = nullptr;
  cf_rel rel = nullptr;
  cf_dget dget = nullptr;
  cf_dcnt dcnt = nullptr;
  cf_dcpy dcpy = nullptr;
  cf_acnt acnt = nullptr;
  cf_aget aget = nullptr;
  cf_tid tid = nullptr;
  cf_atid atid = nullptr;
  io_grp grp = nullptr;
  io_sub sub = nullptr;
  io_smp smp = nullptr;
  io_int ival = nullptr;
  io_str cname = nullptr, unit = nullptr;
  void *subs = nullptr;
  CFMutableDictionaryRef want = nullptr, got = nullptr;
  CFStringRef key = nullptr;
};

// Opened once per process and kept: the subscription is reused by every
// snapshot.
Ior &ior() {
  static Ior r;
  if (r.tried) return r;
  r.tried = true;
  void *cf = dlopen("/System/Library/Frameworks/CoreFoundation.framework/CoreFoundation", RTLD_NOW | RTLD_LOCAL);
  void *io = dlopen("/usr/lib/libIOReport.dylib", RTLD_NOW | RTLD_LOCAL);
  if (cf == nullptr || io == nullptr) {
    r.why = "libIOReport could not be loaded";
    return r;
  }
  r.mks = reinterpret_cast<cf_mks>(dlsym(cf, "CFStringCreateWithCString"));
  r.gets = reinterpret_cast<cf_gets>(dlsym(cf, "CFStringGetCString"));
  r.rel = reinterpret_cast<cf_rel>(dlsym(cf, "CFRelease"));
  r.dget = reinterpret_cast<cf_dget>(dlsym(cf, "CFDictionaryGetValue"));
  r.dcnt = reinterpret_cast<cf_dcnt>(dlsym(cf, "CFDictionaryGetCount"));
  r.dcpy = reinterpret_cast<cf_dcpy>(dlsym(cf, "CFDictionaryCreateMutableCopy"));
  r.acnt = reinterpret_cast<cf_acnt>(dlsym(cf, "CFArrayGetCount"));
  r.aget = reinterpret_cast<cf_aget>(dlsym(cf, "CFArrayGetValueAtIndex"));
  r.tid = reinterpret_cast<cf_tid>(dlsym(cf, "CFGetTypeID"));
  r.atid = reinterpret_cast<cf_atid>(dlsym(cf, "CFArrayGetTypeID"));
  r.grp = reinterpret_cast<io_grp>(dlsym(io, "IOReportCopyChannelsInGroup"));
  r.sub = reinterpret_cast<io_sub>(dlsym(io, "IOReportCreateSubscription"));
  r.smp = reinterpret_cast<io_smp>(dlsym(io, "IOReportCreateSamples"));
  r.ival = reinterpret_cast<io_int>(dlsym(io, "IOReportSimpleGetIntegerValue"));
  r.cname = reinterpret_cast<io_str>(dlsym(io, "IOReportChannelGetChannelName"));
  r.unit = reinterpret_cast<io_str>(dlsym(io, "IOReportChannelGetUnitLabel"));
  if (!r.mks || !r.gets || !r.rel || !r.dget || !r.dcnt || !r.dcpy || !r.acnt || !r.aget || !r.tid || !r.atid ||
      !r.grp || !r.sub || !r.smp || !r.ival || !r.cname || !r.unit) {
    r.why = "libIOReport lacks a function it needs";
    return r;
  }
  CFStringRef g = r.mks(nullptr, "Energy Model", kCFStringEncodingUTF8);
  r.key = r.mks(nullptr, "IOReportChannels", kCFStringEncodingUTF8);
  if (g == nullptr || r.key == nullptr) {
    if (g != nullptr) r.rel(g);
    r.why = "CoreFoundation refused a string";
    return r;
  }
  CFDictionaryRef ch = r.grp(g, nullptr, 0, 0, 0);
  r.rel(g);
  if (ch == nullptr) {
    r.why = "IOReport has no Energy Model channels on this Mac";
    return r;
  }
  r.want = r.dcpy(nullptr, r.dcnt(ch), ch);
  r.rel(ch);
  if (r.want == nullptr) {
    r.why = "CoreFoundation could not copy the channel list";
    return r;
  }
  r.subs = r.sub(nullptr, r.want, &r.got, 0, nullptr);
  if (r.subs == nullptr || r.got == nullptr) {
    r.why = "IOReport refused the subscription";
    return r;
  }
  r.up = true;
  return r;
}

std::string str(Ior &r, CFStringRef s) {
  char b[256] = {0};
  if (s == nullptr || !r.gets(s, b, sizeof b, kCFStringEncodingUTF8)) return "";
  return b;
}

double unit(const std::string &u) {
  if (u == "nJ") return 1e-9;
  if (u == "uJ" || u == "\xc2\xb5J") return 1e-6;
  if (u == "mJ") return 1e-3;
  if (u == "J") return 1;
  return -1;
}

std::string chip() {
  char b[256] = {0};
  std::size_t n = sizeof b - 1;
  if (sysctlbyname("machdep.cpu.brand_string", b, &n, nullptr, 0) != 0) return "Apple silicon";
  return b;
}

// The Energy Model counters into s, or false with the reason.
bool energy(Snap &s, std::string &why) {
  Ior &r = ior();
  if (!r.up) {
    why = r.why;
    return false;
  }
  CFDictionaryRef sm = r.smp(r.subs, r.got, nullptr);
  if (sm == nullptr) {
    why = "IOReport gave no sample";
    return false;
  }
  // channel -> joules; the CPU is "CPU Energy", or per die "DIE_N_CPU Energy"
  std::map<std::string, double> j;
  double cpu = 0;
  int ncpu = 0;
  bool cpu1 = false, badu = false;
  const void *a = r.dget(sm, r.key);
  if (a != nullptr && r.tid(static_cast<CFTypeRef>(a)) == r.atid()) {
    CFArrayRef arr = static_cast<CFArrayRef>(a);
    const CFIndex m = r.acnt(arr);
    for (CFIndex i = 0; i < m; ++i) {
      CFDictionaryRef c = static_cast<CFDictionaryRef>(r.aget(arr, i));
      if (c == nullptr) continue;
      const std::string nm = str(r, r.cname(c));
      const bool die = nm.size() > 11 && nm.compare(nm.size() - 11, 11, "_CPU Energy") == 0;
      if (nm != "CPU Energy" && !die && nm != "GPU Energy" && nm != "ANE" && nm != "DRAM") continue;
      const double f = unit(str(r, r.unit(c)));
      if (f < 0) {
        badu = true;
        continue;
      }
      const double v = static_cast<double>(r.ival(c, 0)) * f;
      if (nm == "CPU Energy") {
        cpu1 = true;
        j[nm] = v;
      } else if (die) {
        cpu += v;
        ++ncpu;
      } else {
        j[nm] = v;
      }
    }
  }
  r.rel(sm);
  if (!cpu1 && ncpu > 0) j["CPU Energy"] = cpu;
  if (j.empty()) {
    why = badu ? "IOReport's energy channels are in a unit not known here" : "IOReport gave no energy channels";
    return false;
  }
  static const char *const nm[] = {"CPU Energy", "GPU Energy", "ANE", "DRAM"};
  static const char *const kd[] = {"cpu", "gpu", "accelerator", "memory"};
  static const char *const id[] = {"ior:cpu", "ior:gpu", "ior:ane", "ior:dram"};
  const std::string ch = chip();
  for (int k = 0; k < 4; ++k) {
    auto it = j.find(nm[k]);
    if (it == j.end()) continue;
    Dev d;
    d.id = id[k];
    d.kind = kd[k];
    d.vendor = "apple";
    d.name = ch + (k == 0 ? " CPU" : k == 1 ? " GPU" : k == 2 ? " neural engine" : " memory");
    d.src = 1;
    d.ctr = it->second;
    d.aon = k == 0 || k == 3;
    s.d.push_back(d);
  }
  return true;
}

std::string qp(const std::string &p) {
  std::string r = "'";
  for (char c : p) r += c == '\'' ? std::string("'\\''") : std::string(1, c);
  return r + "'";
}

}  // namespace

void mac(const Cfg &c, Snap &s) {
  std::string why = "IOReport not used";
  const bool ok = c.on("ioreport") && energy(s, why);
  std::vector<double> gu;
  if (c.on("ioreg") && !c.ioreg.empty()) {
    std::string out;
    int st = 0;
    if (run(qp(c.ioreg) + " -r -d 1 -c IOAccelerator", out, st)) gu = ioacc(out);
  }
  if (ok) {
    for (Dev &d : s.d) {
      if (d.id == "ior:gpu" && !gu.empty()) d.util = gu[0];
    }
  } else if (c.on("ioreport")) {
    Dev d;
    d.id = "ior:cpu";
    d.kind = "cpu";
    d.vendor = "apple";
    d.name = chip() + " CPU";
    d.bad = 1;
    d.why = "not measurable here: " + why;
    s.d.push_back(d);
  }
  if (!ok) {
    for (std::size_t k = 0; k < gu.size(); ++k) {
      Dev d;
      d.id = "ioreg:gpu" + std::to_string(k);
      d.kind = "gpu";
      d.vendor = "apple";
      d.name = "GPU " + std::to_string(k);
      d.util = gu[k];
      d.bad = 1;
      d.why = "power unknown here: utilisation only (" + why + ")";
      s.d.push_back(d);
    }
  }
  if (c.on("battery") && !c.ioreg.empty() && c.on("ioreg")) {
    std::string out, w;
    int st = 0;
    if (run(qp(c.ioreg) + " -rn AppleSmartBattery", out, st)) {
      const double p = iobat(out, w);
      if (w != "no battery found") {
        Dev d;
        d.id = "bat:mac";
        d.kind = "machine";
        d.name = "battery";
        if (p >= 0) {
          d.src = 2;
          d.w = p;
        } else {
          d.bad = 1;
          d.why = w;
        }
        s.d.push_back(d);
      }
    }
  }
  if (c.on("proc") && c.root.empty()) {
    // CPU time in seconds: the host's ticks (user, system, nice; idle) at
    // the clock tick rate, and this process with its reaped children
    host_cpu_load_info_data_t h;
    mach_msg_type_number_t m = HOST_CPU_LOAD_INFO_COUNT;
    const long hz = sysconf(_SC_CLK_TCK);
    static const mach_port_t hp = mach_host_self();
    if (hz > 0 && host_statistics(hp, HOST_CPU_LOAD_INFO, reinterpret_cast<host_info_t>(&h), &m) ==
                      KERN_SUCCESS) {
      const double b = static_cast<double>(h.cpu_ticks[CPU_STATE_USER]) + h.cpu_ticks[CPU_STATE_SYSTEM] +
                       h.cpu_ticks[CPU_STATE_NICE];
      s.scpu = b / hz;
      s.stot = (b + h.cpu_ticks[CPU_STATE_IDLE]) / hz;
    }
    if (c.pid == static_cast<long>(getpid())) {
      rusage a{}, k{};
      if (getrusage(RUSAGE_SELF, &a) == 0 && getrusage(RUSAGE_CHILDREN, &k) == 0) {
        auto sec = [](const timeval &t) { return t.tv_sec + t.tv_usec * 1e-6; };
        s.pcpu = sec(a.ru_utime) + sec(a.ru_stime) + sec(k.ru_utime) + sec(k.ru_stime);
      }
    }
  }
}

#endif

}  // namespace pwr
