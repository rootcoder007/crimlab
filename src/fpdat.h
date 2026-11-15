// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The tables of the footprint family (fpdat.cpp), read by fp.cpp.
#ifndef CLAB_FPDAT_H
#define CLAB_FPDAT_H

namespace cfp {

// One row of the carbon-intensity table: location code, ISO 3166-1
// alpha-3 code, name, lifecycle gCO2e/kWh, year, index into src.
struct Ci {
  const char *loc, *iso3, *name;
  double g;
  int year, src;
};

// An IANA time zone and its territory (ISO 3166-1 alpha-2).
struct Tz {
  const char *tz, *terr;
};

extern const char *const src[];
extern const Ci ci[];
extern const int nci;
extern const Tz tz[];
extern const int ntz;

}  // namespace cfp

#endif
