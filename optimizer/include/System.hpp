#pragma once

#include <cmath>

#include "tagarray.hpp"

#include "SystemData.hpp"

namespace LED {

namespace TA = tagarray;

class System {
public:
  int Norbs;
  double energy;
  std::vector<SystemData> data;
  const std::string name;
  inline System(const std::filesystem::path &filepath) noexcept
      : Norbs(0), energy(0.0), data({}), name(filepath) {
    TA::Container c = TA::Container::load(filepath);
    int64_t nalpha = c["nalpha"]->raw_data<int64_t *>()[0];
    int64_t nbeta = c["nbeta"]->raw_data<int64_t *>()[0];
    this->Norbs = 2 * (nalpha + nbeta);
    int64_t grid_size = c["weight"]->count();
    std::vector<double> rhos_a(grid_size), rhos_b(grid_size), taus_a(grid_size),
        taus_b(grid_size);
    for (int64_t iorb = 0; iorb < nalpha; iorb++) {
      this->energy += c["EIASS"]->raw_data<double *>()[iorb];
      this->energy += c["EIAOS"]->raw_data<double *>()[iorb];
      double *rhoi =
          c["rho.A." + std::to_string(iorb + 1)]->raw_data<double *>();
      double *taui =
          c["tau.A." + std::to_string(iorb + 1)]->raw_data<double *>();
      for (int64_t igrd = 0; igrd < grid_size; igrd++) {
        rhos_a[igrd] += rhoi[igrd];
        taus_a[igrd] += taui[igrd];
      }
    }
    for (int64_t iorb = 0; iorb < nbeta; iorb++) {
      this->energy += c["EIBSS"]->raw_data<double *>()[iorb];
      this->energy += c["EIBOS"]->raw_data<double *>()[iorb];
      double *rhoi =
          c["rho.B." + std::to_string(iorb + 1)]->raw_data<double *>();
      double *taui =
          c["tau.B." + std::to_string(iorb + 1)]->raw_data<double *>();
      for (int64_t igrd = 0; igrd < grid_size; igrd++) {
        rhos_b[igrd] += rhoi[igrd];
        taus_b[igrd] += taui[igrd];
      }
    }
    double *weight = c["weight"]->raw_data<double *>();
    for (int64_t iorb = 0; iorb < nalpha; iorb++) {
      double *rhoi =
          c["rho.A." + std::to_string(iorb + 1)]->raw_data<double *>();
      double *taui =
          c["tau.A." + std::to_string(iorb + 1)]->raw_data<double *>();
      SystemData sd_ass{
          .Npoints = grid_size,
          .energy = c["EIASS"]->raw_data<double *>()[iorb],
          .prefactor = new double[grid_size],
          .features = {new double[grid_size], new double[grid_size],
                       new double[grid_size], new double[grid_size],
                       new double[grid_size]}};
      for (int64_t igrd = 0; igrd < grid_size; igrd++) {
        sd_ass.prefactor[igrd] =
            weight[igrd] * std::sqrt(rhoi[igrd] * (rhos_a[igrd] - rhoi[igrd]));
        sd_ass.features[0][igrd] = rhoi[igrd] / (rhos_a[igrd] - rhoi[igrd]);
        sd_ass.features[1][igrd] = taui[igrd] / (taus_a[igrd] - taui[igrd]);
        sd_ass.features[2][igrd] = taui[igrd] / std::pow(rhoi[igrd], 5. / 3.);
        sd_ass.features[3][igrd] = (taus_a[igrd] - taui[igrd]) /
                                   std::pow(rhos_a[igrd] - rhoi[igrd], 5. / 3.);
        sd_ass.features[4][igrd] = 0.0;
      }
      SystemData sd_aos{
          .Npoints = grid_size,
          .energy = c["EIAOS"]->raw_data<double *>()[iorb],
          .prefactor = new double[grid_size],
          .features = {new double[grid_size], new double[grid_size],
                       new double[grid_size], new double[grid_size],
                       new double[grid_size]}};
      for (int64_t igrd = 0; igrd < grid_size; igrd++) {
        sd_aos.prefactor[igrd] =
            weight[igrd] * std::sqrt(rhoi[igrd] * rhos_b[igrd]);
        sd_aos.features[0][igrd] = rhoi[igrd] / rhos_b[igrd];
        sd_aos.features[1][igrd] = taui[igrd] / taus_b[igrd];
        sd_aos.features[2][igrd] = taui[igrd] / std::pow(rhoi[igrd], 5. / 3.);
        sd_aos.features[3][igrd] =
            taus_b[igrd] / std::pow(rhos_b[igrd], 5. / 3.);
        sd_aos.features[4][igrd] = 0.0;
      }
      data.push_back(sd_ass);
      data.push_back(sd_aos);
    }
    for (int64_t iorb = 0; iorb < nbeta; iorb++) {
      double *rhoi =
          c["rho.B." + std::to_string(iorb + 1)]->raw_data<double *>();
      double *taui =
          c["tau.B." + std::to_string(iorb + 1)]->raw_data<double *>();
      SystemData sd_bss{
          .Npoints = grid_size,
          .energy = c["EIBSS"]->raw_data<double *>()[iorb],
          .prefactor = new double[grid_size],
          .features = {new double[grid_size], new double[grid_size],
                       new double[grid_size], new double[grid_size],
                       new double[grid_size]}};
      for (int64_t igrd = 0; igrd < grid_size; igrd++) {
        sd_bss.prefactor[igrd] =
            weight[igrd] * std::sqrt(rhoi[igrd] * (rhos_b[igrd] - rhoi[igrd]));
        sd_bss.features[0][igrd] = rhoi[igrd] / (rhos_b[igrd] - rhoi[igrd]);
        sd_bss.features[1][igrd] = taui[igrd] / (taus_b[igrd] - taui[igrd]);
        sd_bss.features[2][igrd] = taui[igrd] / std::pow(rhoi[igrd], 5. / 3.);
        sd_bss.features[3][igrd] = (taus_b[igrd] - taui[igrd]) /
                                   std::pow(rhos_b[igrd] - rhoi[igrd], 5. / 3.);
        sd_bss.features[4][igrd] = 0.0;
      }
      SystemData sd_bos{
          .Npoints = grid_size,
          .energy = c["EIBOS"]->raw_data<double *>()[iorb],
          .prefactor = new double[grid_size],
          .features = {new double[grid_size], new double[grid_size],
                       new double[grid_size], new double[grid_size],
                       new double[grid_size]}};
      for (int64_t igrd = 0; igrd < grid_size; igrd++) {
        sd_bos.prefactor[igrd] =
            weight[igrd] * std::sqrt(rhoi[igrd] * rhos_a[igrd]);
        sd_bos.features[0][igrd] = rhoi[igrd] / rhos_b[igrd];
        sd_bos.features[1][igrd] = taui[igrd] / taus_b[igrd];
        sd_bos.features[2][igrd] = taui[igrd] / std::pow(rhoi[igrd], 5. / 3.);
        sd_bos.features[3][igrd] =
            taus_a[igrd] / std::pow(rhos_a[igrd], 5. / 3.);
        sd_bos.features[4][igrd] = 0.0;
      }
      data.push_back(sd_bss);
      data.push_back(sd_bos);
    }
  }
  ~System() {
    for (std::size_t s = 0; s < this->data.size(); s++) {
      delete[] data[s].prefactor;
      for (int i = 0; i < NFEATURES; i++) {
        delete[] data[s].features[i];
      }
    }
  }
};

} // namespace LED
