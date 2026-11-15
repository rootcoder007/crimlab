// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The power and energy readers of hwpwr(): what each device of the machine
// draws, read from the operating system's own counters. A reading is a
// snapshot (Snap) of every device at one instant; the figures over an
// interval come from two or more snapshots (rows()). Each source is read
// through a root folder or a program path the caller can replace, so every
// parser runs on recorded samples without the hardware.
#ifndef CLAB_PWR_H
#define CLAB_PWR_H

#include <string>
#include <vector>

namespace pwr {

// One device at one instant.
struct Dev {
  std::string id;      // stable key within a machine ("rapl:package-0")
  std::string kind;    // cpu, memory, gpu, accelerator, machine, meter
  std::string vendor;  // intel, amd, nvidia, apple, or "" when unknown
  std::string name;
  std::string why;     // why a figure is missing, or a note on it
  int src = 0;         // 0 nothing to read, 1 energy counter, 2 power only
  double ctr = 0;      // the counter, joules (src 1)
  double wrap = 0;     // the counter's range in joules; 0 when it never wraps
  double w = -1;       // power now, watts; < 0 unknown
  double util = -1;    // busy share in [0, 1]; < 0 unknown
  double rated = -1;   // a rated power the machine reports, watts; < 0 none
  double pown = -1;    // the process's share of the device, [0, 1]; < 0 none
  int inpkg = 0;       // 1: its energy is inside the CPU package figure
  int aon = 0;         // 1: an always-on counter (a pause in it is a fault)
  int bad = 0;         // 1: the device is present but gave no number
};

// The whole machine at one instant: t in seconds on a steady clock; the
// process (tree) CPU time pcpu and the system's busy and total CPU time
// (scpu, stot), all in the same unit; < 0 when unknown.
struct Snap {
  double t = 0;
  double pcpu = -1, scpu = -1, stot = -1;
  std::vector<Dev> d;
};

// Where to read, and what to leave out.
struct Cfg {
  std::string root;    // prefix of /sys and /proc ("" the real ones)
  std::string nvml;    // the NVIDIA library ("" its usual names)
  std::string smi;     // nvidia-smi ("" not used)
  std::string ioreg;   // ioreg on macOS ("" not used)
  std::vector<std::string> skip;  // back-ends not used, by name
  long pid = 0;        // the process whose share is reported
  bool on(const char *b) const;
};

// One row of the result.
struct Row {
  std::string kind, vendor, name, scope, how, why;
  double watts, joules, util, flo, fhi;  // NaN when unknown
};

// Rated powers the caller gives (watts; < 0 none), for devices that report a
// utilisation but no power.
struct Rated {
  double cpu = -1, gpu = -1;
};

// A snapshot of every back-end Cfg allows.
Snap snap(const Cfg &c);
// The rows over the snapshots s (time order), with the caller's ratings.
std::vector<Row> rows(const std::vector<Snap> &s, const Rated &r);

// A snapshot as text (one line per device) and back; false on text that is
// not a snapshot.
std::string put(const Snap &s);
bool get(const std::string &txt, Snap &s);

// The back-ends, each adding its devices to s. They are exported for the
// self-test, which drives them with fake trees and recorded text.
void rapl(const Cfg &c, Snap &s);
void drm(const Cfg &c, Snap &s);
void bat(const Cfg &c, Snap &s);
void proc(const Cfg &c, Snap &s);
void nv(const Cfg &c, Snap &s);
void mac(const Cfg &c, Snap &s);
void win(const Cfg &c, Snap &s);

// Parsers of program output. smi(): lines "power, utilisation, name" of
// nvidia-smi --query-gpu=power.draw,utilization.gpu,name
// --format=csv,noheader,nounits. ioacc(): the busy share of each GPU in the
// output of ioreg -r -d 1 -c IOAccelerator. iobat(): the battery's
// discharge in watts from ioreg -rn AppleSmartBattery (< 0 when not on
// battery or not found), with the reason in why.
void smi(const std::string &out, Snap &s);
std::vector<double> ioacc(const std::string &out);
double iobat(const std::string &out, std::string &why);
// The channels of an Energy Meter Interface metadata block (version 1 or
// 2) of n bytes: names, and whether each is in picowatt-hours. False on a
// block that is too short or not of that version.
bool emimeta(const unsigned char *b, unsigned long n, int ver, std::vector<std::string> &names,
             std::vector<int> &pwh);

// Files and programs. rd(): the whole file, false when it cannot be read
// (err set to the system's reason, "denied" for a permission refusal).
bool rd(const std::string &path, std::string &out, std::string &err);
// The output of a command line, false when it could not run or exited
// non-zero (st: its status, -1 when it did not start).
bool run(const std::string &cmd, std::string &out, int &st);
// The entries of a folder, sorted; empty when it cannot be read.
std::vector<std::string> ls(const std::string &dir);

}  // namespace pwr

#endif
