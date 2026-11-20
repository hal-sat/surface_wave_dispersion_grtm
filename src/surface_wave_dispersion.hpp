// SPDX-FileCopyrightText: Copyright 2024, 2026 Harusato Kimura
// SPDX-License-Identifier: MIT

#ifndef surface_wave_dispersion_hpp
#define surface_wave_dispersion_hpp

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

// Velocities and rigidity of every layer at a frequency. They are real for an
// elastic medium and complex for a viscoelastic one.
struct Medium
{
  std::vector<std::complex<double>> vp;
  std::vector<std::complex<double>> vs;
  std::vector<std::complex<double>> mu;

  // List of layer indices that correspond to sandwitched low velocity
  // layers (LVLs), judged by the S-wave velocities of the elastic medium.
  std::vector<int> LVL_indices;
};

class SurfaceWaveDispersion
{
private:
  int nf;
  double fmin;
  double fmax;
  double df;
  // freq[i] = fmin + df * i for i = 0 to nf, where df = (fmax - fmin) / nf.

  int ngrid = 500;
  // Number of grid points for grid search before bisection is constructed by
  // ngrid * predicted_number_of_modes(), where predicted_number_of_modes() is
  // given by Fan et al. (2007, CJG).

  std::vector<double> freq;
  std::vector<std::vector<double>> phvs_rayleigh;
  // keys_rayleigh[i][j] / 2 is the index of the secular function with which
  // phvs_rayleigh[i][j] was found (0: free surface, 1, 2, ...: LVLs).
  std::vector<std::vector<int>> keys_rayleigh;
  std::vector<std::vector<double>> medium_responses_rayleigh;
  // medium_responses_rayleigh[i][j]: absolute value of the medium response
  // of the jth mode at freq[i] (vertical point force and vertical
  // displacement at the free surface).
  std::vector<std::vector<double>> psds_rayleigh;
  // psds_rayleigh[i][j]: power of the jth mode relative to that of the 0th
  // mode at freq[i] (vertical point force and vertical displacement at the
  // free surface).
  std::vector<std::vector<double>> attens_rayleigh;
  // attens_rayleigh[i][j]: attenuation coefficient alpha (1/m) of the jth
  // mode at freq[i]; the wavenumber of the mode is omega / phv - i * alpha.
  size_t nmodes_rayleigh = 1;
  double fluctuation_for_root_search = 1e-5;

  std::vector<double> thickness;
  std::vector<double> rho;

  // Vp and Vs of the profile and, if given, Qp and Qs. Q <= 0 means that the
  // layer is elastic.
  std::vector<double> vp_ref;
  std::vector<double> vs_ref;
  std::vector<double> qp;
  std::vector<double> qs;
  bool attenuation_given = false;
  // If true, Qp and Qs of the profile are ignored and the medium is treated
  // as elastic.
  bool attenuation_ignored = false;

  // Reference frequency (Hz) at which Vp and Vs of the profile are defined.
  // If it is given, the causal constant-Q model of Kjartansson (1979) is
  // used; otherwise the complex velocities are frequency independent.
  bool reference_frequency_given = false;
  double f_ref = 0.0;
  // Exponents of the power law of the constant-Q model.
  std::vector<double> gamma_p;
  std::vector<double> gamma_s;

  // Elastic and viscoelastic media at the reference frequency (at all the
  // frequencies if the reference frequency is not given).
  Medium medium_elastic_ref;
  Medium medium_viscoelastic_ref;

  std::vector<double> depth;
  // depth[0] = z_{0}
  // Caution!!! depth[i] = z_{ind} but
  // var[i] = var_{ind - 1} for var in thickness, vp, vs, rho, mu
  // thickness[j] = depth[j + 1] - depth[j]
  // for (j = ind - 1)th layer, the depth of bottom interface is depth[j + 1]
  // and that of above interface is depth[j].

  int const bisect_max_itr = 100;
  double const bisect_acc = 1e-3;
  // A sign change is regarded as not being a root (a pole or a jump of the
  // secular function) if the secular function at the middle of the interval
  // is larger than bisect_reject_factor times what a root allows; see
  // bisection_sf().
  double const bisect_reject_factor = 100.0;

  // Parameters of the contour integration around each root in the complex
  // wavenumber plane (used for the medium response). The radius of the
  // circle is the smaller of contour_rel_radius * k and
  // contour_gap_fraction * (distance to the nearest other root).
  int const contour_npts = 32;
  double const contour_rel_radius = 2e-4;
  double const contour_gap_fraction = 0.25;

  // Parameters of the search for the complex roots of a viscoelastic medium.
  int const newton_max_itr = 50;
  double const newton_rel_tol = 1e-13;
  // Upper limit of a Newton step relative to the wavenumber.
  double const newton_rel_step_limit = 0.02;
  // Step of the finite difference relative to the wavenumber.
  double const newton_rel_diff_step = 1e-7;
  // A root is accepted if its real part is within
  // (1 +- wavenumber_rel_range) * (wavenumber of the elastic medium) and if
  // its attenuation coefficient is not negative and below alpha_max (1/m).
  double const wavenumber_rel_range = 0.2;
  double const alpha_max = 1.0;
  // Two roots are regarded as the same one if their wavenumbers agree
  // within this relative tolerance.
  double const coincidence_rel_tol = 1e-6;
  // Number of steps in which the attenuation is switched on to follow a
  // root from the elastic medium.
  int const tracking_num_steps = 10;

  void build_reference_media();
  void thickness2depth();
  void inspect_nmodes_rayleigh();
  void generate_frequency_vector();

  // Medium at the angular frequency omega. Without the reference frequency
  // it is the one at the reference frequency for any omega.
  Medium medium_at(double const omega, bool const viscoelastic);

  // Medium between the elastic one (fraction = 0) and the viscoelastic one
  // (fraction = 1); the complex velocities are interpolated geometrically.
  Medium interpolate_media(Medium const& elastic,
                           Medium const& viscoelastic,
                           double const fraction);

  // returns F_0, F_l1, F_l2, ... for a (possibly complex) horizontal
  // wavenumber k, where l1, l2, ... are the layer indices that correspond to
  // LVLs. F_0 is the secular function of the free surface and F_l is that of
  // the LVL l.
  std::vector<std::complex<double>> rayleigh_secular_functions(
    std::complex<double> const k,
    double const omega,
    Medium const& medium);

  // returns Re[F_0], Im[F_0], Re[F_l1], Im[F_l1], Re[F_l2], Im[F_l2], ...
  // for the real wavenumber omega / phv.
  std::vector<double> rayleigh_secular_function(double const phv,
                                                double const omega,
                                                Medium const& medium);

  // Prepares E, E^{-1} and the Lambda kernel of every layer and computes the
  // generalized reflection coefficients \hat{R}_{du}^{(j)} (j = N + 1 to 1)
  // for a (possibly complex) horizontal wavenumber k.
  // vec_R_du_i[i] = \hat{R}_{du}^{(i + 1)}.
  void rayleigh_downward_recursion(
    std::complex<double> const k,
    double const omega,
    Medium const& medium,
    std::vector<Eigen::Matrix<std::complex<double>, 4, 4>>& vec_E_i,
    std::vector<Eigen::Matrix<std::complex<double>, 4, 4>>& vec_E_inv_i,
    std::vector<Eigen::Matrix<std::complex<double>, 2, 2>>& vec_Lambda_kernel_i,
    std::vector<Eigen::Matrix<std::complex<double>, 2, 2>>& vec_R_du_i);

  // Surface response Y(omega, k): (displacement at the free surface) =
  // Y * (traction at the free surface). Rows and columns are ordered as
  // (horizontal, vertical). Its poles are the Rayleigh-wave modes.
  Eigen::Matrix<std::complex<double>, 2, 2> rayleigh_surface_response(
    std::complex<double> const k,
    double const omega,
    Medium const& medium);

  // Medium response A = k_n * Res_{k = k_n} Y_zz(k) of the mode whose root is
  // k_n = omega / phv, evaluated by a contour integral of given radius.
  std::complex<double> rayleigh_medium_response(double const phv,
                                                double const omega,
                                                double const radius,
                                                Medium const& medium);

  void write_table(std::string const file_path,
                   std::vector<std::vector<double>> const& table,
                   size_t const nmodes);

  // v_min and v_max are rayleigh_secular_function() at x_min and x_max.
  double bisection_sf(double const x_min,
                      double const x_max,
                      std::vector<double> const& v_min,
                      std::vector<double> const& v_max,
                      double const omega,
                      int const key,
                      Medium const& medium);

  double vector_rms(std::vector<double> const& vec);
  double vector_rms(std::vector<double> const& vec, size_t ind_s, size_t ind_e);

  double predicted_num_of_modes(double phv, double omega, Medium const& medium);

  bool equals_vs_of_a_layer(double phv, Medium const& medium);
  bool equals_vp_of_a_layer(double phv, Medium const& medium);

  std::vector<int> detect_LVLs(std::vector<std::complex<double>> const& vs);
  bool is_LVL(int const ind, Medium const& medium);

  // --- Complex roots of a viscoelastic medium

  // Newton iteration in the complex wavenumber plane for a root of the
  // sf_index-th secular function, starting from k. Returns false if it does
  // not converge.
  bool newton_root(std::complex<double>& k,
                   double const omega,
                   int const sf_index,
                   Medium const& medium);

  // Whether the complex wavenumber k is acceptable as the root that
  // originates from the wavenumber k_elastic of the elastic medium.
  bool is_acceptable_root(std::complex<double> const k, double const k_elastic);

  // Initial guess of a root of the viscoelastic medium, obtained by
  // following the root k of the elastic medium while the attenuation is
  // switched on step by step.
  std::complex<double> track_root_from_elastic(std::complex<double> k,
                                               double const omega,
                                               int const sf_index,
                                               Medium const& elastic,
                                               Medium const& viscoelastic);

  // Fills the frequencies at which a mode is missing after its cut-in by a
  // linear interpolation along the frequency. table[i][j] is the phase
  // velocity of the jth mode at freq[i] (nan if it does not exist).
  void interpolate_missing_modes(std::vector<std::vector<double>>& table);

  // ---
  // For the sake of simplicity in the computation of secular function.
  std::complex<double> nu(std::complex<double> const k,
                          double const omega,
                          int const ind,
                          Medium const& medium);
  std::complex<double> gamma(std::complex<double> const k,
                             double const omega,
                             int const ind,
                             Medium const& medium);
  std::complex<double> chi1(std::complex<double> const k,
                            double const omega,
                            int const ind,
                            Medium const& medium);
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
                        double const new_fmin,
                        double const new_fmax);

  size_t get_number_of_modes_rayleigh();

  // Whether Qp and Qs are given in the profile (and not ignored).
  bool has_attenuation();

  // Treat the medium as elastic even if Qp and Qs are given in the profile.
  // This must be called before read_profile().
  void set_ignore_attenuation(bool const new_attenuation_ignored);

  void set_number_of_grids(int const new_ngrid);

  void set_basic_parameters(int const new_nf,
                            double const new_fmin,
                            double const new_fmax);

  // Use the constant-Q model of Kjartansson (1979) with the velocities of
  // the profile defined at new_f_ref (Hz).
  void set_reference_frequency(double const new_f_ref);

  void dump_basic_parameters();
  void dump_ground_structure();
  void dump_number_of_modes_rayleigh();
  void dump_LVL_indices();

  void write_rayleigh_results(std::string const file_path);
  void write_rayleigh_power_fraction(std::string const file_path);
  void write_rayleigh_medium_response(std::string const file_path);
  void write_rayleigh_attenuation(std::string const file_path);

  void read_profile(std::string const path, char const delimiter = ',');
  // Phase velocities of the elastic medium.
  void calculate_rayleigh_phase_velocities();
  void calculate_rayleigh_power_fraction();
  // Phase velocities and attenuation coefficients of the viscoelastic
  // medium, starting from the phase velocities of the elastic medium.
  // calculate_rayleigh_phase_velocities() must be called in advance.
  void calculate_rayleigh_attenuation();
};

#endif
