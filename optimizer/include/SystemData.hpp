#pragma once

#include <cstdint>

#include "Config.hpp"

namespace LED {

struct SystemData {
  int64_t Npoints;
  double energy;
  double *prefactor;
  double *features[NFEATURES];
};

} // namespace LED
