#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Utility.h"
#include <map>

/// This namespace holds some useful constants
namespace constants
{
  /// T0 calibration is taken from lab measurement
  const std::vector<std::vector<Float_t>> dirichTimeCorr {
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {}
  };

  const std::map<DiRICH, std::string> dirichNameMap {
    { DiRICH::d1200, "1200" },
    { DiRICH::d1201, "1201" },
    { DiRICH::d1202, "1202" },
    { DiRICH::d1203, "1203" },
    { DiRICH::d1204, "1204" },
    { DiRICH::d1205, "1205" },
    { DiRICH::d1206, "1206" },
    { DiRICH::d1207, "1207" },
    { DiRICH::d1208, "1208" },
    { DiRICH::d1209, "1209" },
    { DiRICH::d120a, "120a" },
    { DiRICH::d120b, "120b" }
  };

  const std::vector<std::string> dirichNames{ "1200",
                                              "1201",
                                              "1202",
                                              "1203",
                                              "1204",
                                              "1205",
                                              "1206",
                                              "1207",
                                              "1208",
                                              "1209",
                                              "120a",
                                              "120b"
  };

} // namespace constants

#endif