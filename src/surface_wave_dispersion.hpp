// SPDX-FileCopyrightText: Copyright 2024 Harusato Kimura

#ifndef surface_wave_dispersion_hpp
#define surface_wave_dispersion_hpp

#include <Eigen/Core>
#include <Eigen/Dense>

#include <cerrno>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

class SurfaceWaveDispersion
{
private:
  int nf;
  double fmax;
  double df;

  int ngrid = 500;
  // Number of grid points for grid search before bisection is constructed by
  // ngrid * predicted_number_of_modes(), where predicted_number_of_modes() is
  // given by Fan et al. (2007, CJG).

  std::vector<double> freq;
  std::vector<std::vector<double>> phvs_love;
  std::vector<std::vector<double>> phvs_rayleigh;
  size_t nmodes_love = 1;
  size_t nmodes_rayleigh = 1;
  double fluctuation_for_root_search = 1e-5;

  std::vector<double> thickness;
  std::vector<double> vp;
  std::vector<double> vs;
  std::vector<double> rho;
  std::vector<double> mu;

  std::vector<double> depth;
  // depth[0] = z_{0}
  // Caution!!! depth[i] = z_{ind} but
  // var[i] = var_{ind - 1} for var in thickness, vp, vs, rho, mu
  // thickness[j] = depth[j + 1] - depth[j]
  // for (j = ind - 1)th layer, the depth of bottom interface is depth[j + 1]
  // and that of above interface is depth[j].

  int const bisect_max_itr = 100;
  double const bisect_acc = 1e-3;

  std::vector<int> LVL_indices;
  // List of layer indices that correspond to sandwitched low velocity
  // layers (LVLs).

  void calc_mu();
  void thickness2depth();
  void inspect_nmodes_love();
  void inspect_nmodes_rayleigh();
  void generate_frequency_vector();

  // returns Re[F_0], Im[F_0], Re[F_l1], Im[F_l1], Re[F_l2], Im[F_l2], ...
  // where l1, l2, ... are the layer indices that correspond to LVLs.
  // std::complex<double> rayleigh_secular_function(double const phv,
  //                                                double const omega);
  std::vector<double> rayleigh_secular_function(double const phv,
                                                double const omega);
  std::complex<double> love_secular_function(double const phv,
                                             double const omega);

  double bisection_sf(double const x_min,
                      double const x_max,
                      double const omega,
                      int const key,
                      std::string target);

  double vector_rms(std::vector<double> const& vec);
  double vector_rms(std::vector<double> const& vec, size_t ind_s, size_t ind_e);

  double predicted_num_of_modes(double phv, double omega);

  bool equals_vs_of_a_layer(double phv);
  bool equals_vp_of_a_layer(double phv);

  void detect_LVLs();
  bool is_LVL(int const ind);

  // ---
  // For the sake of simplicity in the computation of secular function.
  std::complex<double> nu(std::complex<double> const k,
                          double const omega,
                          int const ind);
  std::complex<double> gamma(std::complex<double> const k,
                             double const omega,
                             int const ind);
  std::complex<double> chi1(std::complex<double> const k,
                            double const omega,
                            int const ind);
  Eigen::Matrix<std::complex<double>, 2, 2> Lambda_d(
    int const ind,
    double const z,
    std::complex<double> const gamma_ind,
    std::complex<double> const nu_ind);
  Eigen::Matrix<std::complex<double>, 2, 2> Lambda_u(
    int const ind,
    double const z,
    std::complex<double> const gamma_ind,
    std::complex<double> const nu_ind);

  // ---

public:
  SurfaceWaveDispersion(int const new_nf,
                        double const new_fmax,
                        double const new_df);

  size_t get_number_of_modes_love();
  size_t get_number_of_modes_rayleigh();

  void set_number_of_grids(int const new_ngrid);

  void set_basic_parameters(int const new_nf,
                            double const new_fmax,
                            double const new_df);

  void dump_basic_parameters();
  void dump_ground_structure();
  void dump_number_of_modes_love();
  void dump_number_of_modes_rayleigh();
  void dump_LVL_indices();

  void write_love_results(std::string const file_path);
  void write_rayleigh_results(std::string const file_path);

  void read_profile(std::string const path, char const delimiter = ',');
  void calculate_love_phase_velocities();
  void calculate_rayleigh_phase_velocities();
  void calculate_love_power_fraction();
  void calculate_rayleigh_power_fraction();
};

#endif
