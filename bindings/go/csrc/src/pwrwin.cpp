// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The Windows readers of hwpwr() (src/pwr.h), on documented interfaces
// only. The CPU packages: the Energy Meter Interface (EMI) devices, found
// by their device interface class and read with IOCTL_EMI_GET_VERSION,
// _GET_METADATA_SIZE, _GET_METADATA and _GET_MEASUREMENT (absolute energy
// in picowatt-hours). Channels whose name says PKG are CPU packages and
// those that say DRAM memory; any other channel is listed and not added.
// Whether a machine has an EMI device, and whether an ordinary user may
// open it, varies: each case is reported as found. The battery: its
// discharge rate in milliwatts from IOCTL_BATTERY_QUERY_STATUS. The CPU
// time of this process and of the system: GetProcessTimes() and
// GetSystemTimes(). The configuration manager (cfgmgr32.dll) is loaded at
// run time, so nothing beyond the system's own kernel32 is linked.
#include "pwr.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <winioctl.h>
#include <cfgmgr32.h>
// emi.h and batclass.h are in the Windows SDK; the MinGW headers that R
// is built with lack the first and may lack the second. Where a header is
// missing, the few names used here are declared as the SDK documents
// them (Energy Meter Interface; IOCTL_BATTERY_QUERY_STATUS).
#if defined(__has_include) && __has_include(<emi.h>)
#include <emi.h>
#else
#define EMI_VERSION_V1 1
#define EMI_VERSION_V2 2
#define IOCTL_EMI_GET_VERSION CTL_CODE(FILE_DEVICE_UNKNOWN, 0x0, METHOD_BUFFERED, FILE_READ_ACCESS)
#define IOCTL_EMI_GET_METADATA_SIZE CTL_CODE(FILE_DEVICE_UNKNOWN, 0x1, METHOD_BUFFERED, FILE_READ_ACCESS)
#define IOCTL_EMI_GET_METADATA CTL_CODE(FILE_DEVICE_UNKNOWN, 0x2, METHOD_BUFFERED, FILE_READ_ACCESS)
#define IOCTL_EMI_GET_MEASUREMENT CTL_CODE(FILE_DEVICE_UNKNOWN, 0x3, METHOD_BUFFERED, FILE_READ_ACCESS)
struct EMI_VERSION {
  USHORT EmiVersion;
};
struct EMI_METADATA_SIZE {
  ULONG MetadataSize;
};
struct EMI_CHANNEL_MEASUREMENT_DATA {
  ULONGLONG AbsoluteEnergy;
  ULONGLONG AbsoluteTime;
};
#endif
#if defined(__has_include) && __has_include(<batclass.h>)
#include <batclass.h>
#else
#ifndef FILE_DEVICE_BATTERY
#define FILE_DEVICE_BATTERY 0x00000029
#endif
#define IOCTL_BATTERY_QUERY_TAG CTL_CODE(FILE_DEVICE_BATTERY, 0x10, METHOD_BUFFERED, FILE_READ_ACCESS)
#define IOCTL_BATTERY_QUERY_STATUS CTL_CODE(FILE_DEVICE_BATTERY, 0x13, METHOD_BUFFERED, FILE_READ_ACCESS)
#define BATTERY_DISCHARGING 0x00000002
#define BATTERY_UNKNOWN_RATE 0x80000000
struct BATTERY_WAIT_STATUS {
  ULONG BatteryTag;
  ULONG Timeout;
  ULONG PowerState;
  ULONG LowCapacity;
  ULONG HighCapacity;
};
struct BATTERY_STATUS {
  ULONG PowerState;
  ULONG Capacity;
  ULONG Voltage;
  LONG Rate;
};
#endif

#include <cmath>
#include <cstring>
#include <cwchar>
#include <string>
#include <vector>
#endif

namespace pwr {

#ifndef _WIN32

void win(const Cfg &, Snap &) {}

#else

namespace {

// {45BD8344-7ED6-49CF-A440-C276C933B053}, GUID_DEVICE_ENERGY_METER
const GUID gemi = {0x45bd8344, 0x7ed6, 0x49cf, {0xa4, 0x40, 0xc2, 0x76, 0xc9, 0x33, 0xb0, 0x53}};
// {72631E54-78A4-11D0-BCF7-00AA00B7B32A}, GUID_DEVCLASS_BATTERY
const GUID gbat = {0x72631e54, 0x78a4, 0x11d0, {0xbc, 0xf7, 0x00, 0xaa, 0x00, 0xb7, 0xb3, 0x2a}};

typedef CONFIGRET(WINAPI *cm_sz)(PULONG, LPGUID, DEVINSTID_W, ULONG);
typedef CONFIGRET(WINAPI *cm_ls)(LPGUID, DEVINSTID_W, PZZWSTR, ULONG, ULONG);

// The paths of the present interfaces of class g.
std::vector<std::wstring> ifaces(const GUID &g) {
  std::vector<std::wstring> v;
  static HMODULE cm = LoadLibraryW(L"cfgmgr32.dll");
  if (cm == nullptr) return v;
  cm_sz sz = reinterpret_cast<cm_sz>(reinterpret_cast<void *>(GetProcAddress(cm, "CM_Get_Device_Interface_List_SizeW")));
  cm_ls ls = reinterpret_cast<cm_ls>(reinterpret_cast<void *>(GetProcAddress(cm, "CM_Get_Device_Interface_ListW")));
  if (sz == nullptr || ls == nullptr) return v;
  GUID gg = g;
  for (int tries = 0; tries < 3; ++tries) {
    ULONG n = 0;
    if (sz(&n, &gg, nullptr, CM_GET_DEVICE_INTERFACE_LIST_PRESENT) != CR_SUCCESS || n == 0) return v;
    std::vector<wchar_t> b(n + 1, 0);
    const CONFIGRET r = ls(&gg, nullptr, b.data(), n, CM_GET_DEVICE_INTERFACE_LIST_PRESENT);
    if (r == CR_BUFFER_SMALL) continue;
    if (r != CR_SUCCESS) return v;
    for (const wchar_t *p = b.data(); *p != 0; p += std::wcslen(p) + 1) v.emplace_back(p);
    return v;
  }
  return v;
}

bool ioc(HANDLE h, DWORD code, void *in, DWORD nin, void *out, DWORD nout, DWORD &got) {
  got = 0;
  return DeviceIoControl(h, code, in, nin, out, nout, &got, nullptr) != 0;
}

std::string s8(const std::wstring &w) {
  std::string s;
  for (wchar_t c : w) s += c < 128 ? static_cast<char>(c) : '?';
  return s;
}

void emi(Snap &s) {
  const std::vector<std::wstring> dv = ifaces(gemi);
  bool any = false;
  int k = 0;
  std::string why;
  for (const std::wstring &p : dv) {
    ++k;
    HANDLE h = CreateFileW(p.c_str(), GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, OPEN_EXISTING, 0,
                           nullptr);
    if (h == INVALID_HANDLE_VALUE) {
      why = GetLastError() == ERROR_ACCESS_DENIED ? "an energy meter (EMI) is present but not readable: access denied"
                                                  : "an energy meter (EMI) is present but could not be opened";
      continue;
    }
    DWORD got = 0;
    EMI_VERSION ver{};
    EMI_METADATA_SIZE msz{};
    std::vector<unsigned char> meta;
    std::vector<std::string> nm;
    std::vector<int> pwh;
    bool ok = ioc(h, IOCTL_EMI_GET_VERSION, nullptr, 0, &ver, sizeof ver, got) && got >= sizeof ver &&
              (ver.EmiVersion == EMI_VERSION_V1 || ver.EmiVersion == EMI_VERSION_V2) &&
              ioc(h, IOCTL_EMI_GET_METADATA_SIZE, nullptr, 0, &msz, sizeof msz, got) && got >= sizeof msz &&
              msz.MetadataSize > 0 && msz.MetadataSize < (1u << 20);
    if (ok) {
      meta.assign(msz.MetadataSize, 0);
      ok = ioc(h, IOCTL_EMI_GET_METADATA, nullptr, 0, meta.data(), msz.MetadataSize, got) &&
           emimeta(meta.data(), got, ver.EmiVersion, nm, pwh) && !nm.empty();
    }
    std::vector<EMI_CHANNEL_MEASUREMENT_DATA> m(nm.empty() ? 1 : nm.size());
    if (ok) {
      const DWORD want = static_cast<DWORD>(m.size() * sizeof(EMI_CHANNEL_MEASUREMENT_DATA));
      ok = ioc(h, IOCTL_EMI_GET_MEASUREMENT, nullptr, 0, m.data(), want, got) && got >= want;
    }
    CloseHandle(h);
    if (!ok) {
      why = "an energy meter (EMI) is present but did not answer as documented";
      continue;
    }
    any = true;
    for (std::size_t i = 0; i < nm.size(); ++i) {
      Dev d;
      d.id = "emi:" + std::to_string(k) + ":" + nm[i];
      const bool pkg = nm[i].find("PKG") != std::string::npos, mem = nm[i].find("DRAM") != std::string::npos;
      d.kind = pkg ? "cpu" : mem ? "memory" : "meter";
      d.name = nm[i];
      if (!pwh[i]) {
        d.bad = 1;
        d.why = "the meter's unit is not picowatt-hours";
      } else {
        d.src = 1;
        // 1 pWh = 3.6e-9 J
        d.ctr = static_cast<double>(m[i].AbsoluteEnergy) * 3.6e-9;
      }
      if (!pkg && !mem) d.why = "a meter channel of unknown scope: listed, not added";
      s.d.push_back(d);
    }
  }
  if (!any) {
    Dev d;
    d.id = "emi:none";
    d.kind = "cpu";
    d.name = "CPU";
    d.bad = 1;
    d.why = !why.empty() ? why
                         : "not measurable here: no energy meter device (EMI); Windows has one only on some "
                           "machines and not in virtual machines";
    s.d.push_back(d);
  }
}

void battery(Snap &s) {
  int k = 0;
  for (const std::wstring &p : ifaces(gbat)) {
    HANDLE h = CreateFileW(p.c_str(), GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr,
                           OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (h == INVALID_HANDLE_VALUE) continue;
    DWORD got = 0;
    ULONG wait = 0, tag = 0;
    BATTERY_WAIT_STATUS ws{};
    BATTERY_STATUS bs{};
    const bool ok = ioc(h, IOCTL_BATTERY_QUERY_TAG, &wait, sizeof wait, &tag, sizeof tag, got) && tag != 0 &&
                    ((ws.BatteryTag = tag), true) &&
                    ioc(h, IOCTL_BATTERY_QUERY_STATUS, &ws, sizeof ws, &bs, sizeof bs, got) && got >= sizeof bs;
    CloseHandle(h);
    if (!ok) continue;
    Dev d;
    d.id = "bat:" + std::to_string(k++);
    d.kind = "machine";
    d.name = "battery " + std::to_string(k);
    if (!(bs.PowerState & BATTERY_DISCHARGING)) {
      d.bad = 1;
      d.why = "not on battery: the battery gives the machine's power only while it discharges";
    } else if (static_cast<ULONG>(bs.Rate) == BATTERY_UNKNOWN_RATE || bs.Rate == 0) {
      d.bad = 1;
      d.why = "the battery reports no discharge rate";
    } else {
      d.src = 2;
      d.w = std::fabs(static_cast<double>(bs.Rate)) * 1e-3;
    }
    s.d.push_back(d);
  }
}

double ft(const FILETIME &f) {
  return (static_cast<double>(f.dwHighDateTime) * 4294967296.0 + f.dwLowDateTime) * 1e-7;
}

}  // namespace

void win(const Cfg &c, Snap &s) {
  if (c.on("emi")) emi(s);
  if (c.on("battery")) battery(s);
  if (c.on("proc") && c.root.empty()) {
    FILETIME id, kn, us, cr, ex, pk, pu;
    // kernel time includes idle time
    if (GetSystemTimes(&id, &kn, &us)) {
      s.scpu = ft(kn) - ft(id) + ft(us);
      s.stot = ft(kn) + ft(us);
    }
    if (c.pid == static_cast<long>(GetCurrentProcessId()) && GetProcessTimes(GetCurrentProcess(), &cr, &ex, &pk, &pu)) {
      s.pcpu = ft(pk) + ft(pu);
    }
  }
}

#endif

}  // namespace pwr
