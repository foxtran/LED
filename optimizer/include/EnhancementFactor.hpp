#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <tuple>
#include <unordered_map>

#include "System.hpp"

namespace LED {

using ComputeError =
    std::function<double(const System &, const double *const, const double)>;
using KeyComputeError = std::tuple<int, int, bool, bool>;

struct KeyComputeError_hash {
  std::size_t operator()(const KeyComputeError &k) const {
    return (std::get<0>(k) << 24) + (std::get<1>(k) << 16) +
           ((std::size_t)std::get<2>(k) << 8) +
           ((std::size_t)std::get<3>(k) << 7);
  }
};

constexpr int get_features_count(bool use_rho, bool use_tau, bool use_nrg) {
  const bool rho_used[NFEATURES] = {true, false, true, true, true};
  const bool tau_used[NFEATURES] = {false, true, true, true, false};
  const bool nrg_used[NFEATURES] = {false, false, false, false, true};
  int vars_count = 0;
  for (int i = 0; i < NFEATURES; i++) {
    if ((use_rho || !rho_used[i]) && (use_tau || !tau_used[i]) &&
        (use_nrg || !nrg_used[i])) {
      vars_count++;
    }
  }
  return vars_count;
}

template <std::size_t GFlevel, bool use_tau = false, bool use_nrg = false>
inline double gorner_factorization_evaluate(const SystemData &data,
                                            const int datapoint,
                                            const double *const coeffs) {
  constexpr bool use_rho = true;
  constexpr bool rho_used[NFEATURES] = {true, false, true, true, true};
  constexpr bool tau_used[NFEATURES] = {false, true, true, true, false};
  constexpr bool nrg_used[NFEATURES] = {false, false, false, false, true};
  constexpr int vars_count = get_features_count(use_rho, use_tau, use_nrg);
  std::array<double, vars_count> v, f;
  for (int i = 0, j = 0; i < NFEATURES; i++) {
    if ((use_rho || !rho_used[i]) && (use_tau || !tau_used[i]) &&
        (use_nrg || !nrg_used[i])) {
      v[j] = data.features[i][datapoint];
      j++;
    }
  }
  for (int i = 0; i < vars_count; i++) {
    f[i] = 1.0;
    for (int j = GFlevel; j >= 1; j--) {
      f[i] = f[i] * coeffs[j + i * GFlevel] * v[i] + 1.;
    }
  }
  double enhancement_factor = 1.0;
  for (int i = 0; i < vars_count; i++) {
    enhancement_factor *= f[i];
  }
  return enhancement_factor;
}

template <double (*compute)(const SystemData &, const int, const double *const)>
double compute_orbital_energy(const SystemData &data,
                              const double *const coeffs) {
  double energy = 0.0;
  for (int i = 0; i < data.Npoints; i++) {
    if (std::abs(data.prefactor[i]) > 1e-15) {
      energy += data.prefactor[i] * coeffs[0] * compute(data, i, &coeffs[0]);
    }
  }
  return energy;
}

template <double (*compute_orbital_energy)(const SystemData &data,
                                           const double *const coeffs)>
double compute_system_error(const System &sys, const double *const coeffs,
                            const double omega) {
  double error = 0.0;
  double total_energy = 0.0;
  for (int i = 0; i < sys.Norbs; i++) {
    double energy = compute_orbital_energy(sys.data[i], coeffs);
    double diff = std::abs(sys.data[i].energy - energy);
    error += omega * diff * diff;
    total_energy += energy;
  }
  error += (1. - omega) * std::abs(total_energy - sys.energy);
  return error;
}

template <ssize_t N>
void fill_map(std::unordered_map<KeyComputeError, ComputeError,
                                 KeyComputeError_hash> &ce) {
  if constexpr (N >= 0) {
    fill_map<N - 1>(ce);
    ce[{1, N, false, false}] = compute_system_error<
        compute_orbital_energy<gorner_factorization_evaluate<N, false, false>>>;
    ce[{1, N, false, true}] = compute_system_error<
        compute_orbital_energy<gorner_factorization_evaluate<N, false, true>>>;
    ce[{1, N, true, false}] = compute_system_error<
        compute_orbital_energy<gorner_factorization_evaluate<N, true, false>>>;
    ce[{1, N, true, true}] = compute_system_error<
        compute_orbital_energy<gorner_factorization_evaluate<N, true, true>>>;
  }
}

class EnhancementFactor {
private:
  ComputeError compute_error_f;
  std::vector<double> coeffs;
  std::size_t Nparams;

public:
  EnhancementFactor(int factorization_order, bool use_tau, bool use_nrg) {
    std::unordered_map<KeyComputeError, ComputeError, KeyComputeError_hash>
        ce_funs = {};
    fill_map<MAX_FACTORIZATION>(ce_funs);
    if (factorization_order > MAX_FACTORIZATION) {
      throw std::runtime_error("Factorization order is higher than " +
                               std::to_string(MAX_FACTORIZATION) + "!");
    }
    if (factorization_order < 0) {
      throw std::runtime_error("Factorization order is less than 0!");
    }
    this->compute_error_f = ce_funs[{1, factorization_order, use_tau, use_nrg}];
    this->Nparams =
        1 + get_features_count(true, use_tau, use_nrg) * factorization_order;
    coeffs = std::vector<double>(this->Nparams, 0.0);
  };
  void set_coefficients(const std::vector<double> &coeffs) {
    if (coeffs.size() != this->Nparams) {
      throw std::runtime_error("Insufficient number of coefficients: " +
                               std::to_string(coeffs.size()) +
                               ". It must be: " + std::to_string(Nparams) +
                               "!");
    }
    this->coeffs = coeffs;
  }
  const std::vector<double> &get_coefficients() noexcept {
    return this->coeffs;
  }
  int get_N_free_parameters() noexcept { return this->Nparams; }
  double compute_error(const System &sys, double omega = 0.5) noexcept {
    return this->compute_error_f(sys, this->coeffs.data(), omega);
  }
  double compute_error(const std::vector<System> &syss,
                       double omega = 0.5) noexcept {
    double res = 0.0;
    for (std::size_t s = 0; s < syss.size(); s++) {
      res += this->compute_error_f(syss[s], this->coeffs.data(), omega);
    }
    return res;
  }
};

} // namespace LED
