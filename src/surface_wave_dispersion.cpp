// SPDX-FileCopyrightText: Copyright 2024, 2026 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "surface_wave_dispersion.hpp"

SurfaceWaveDispersion::SurfaceWaveDispersion(int const new_nf,
                                             double const new_fmax,
                                             double const new_df)
{
  this->nf = new_nf;
  this->fmax = new_fmax;
  this->df = new_df;

  this->generate_frequency_vector();
}

void
SurfaceWaveDispersion::generate_frequency_vector()
{
  this->freq.clear();
  this->freq.resize(this->nf + 1);
  for (int i = 0; i < this->nf + 1; ++i) {
    this->freq[i] = this->fmax / (double)this->nf * (double)i;
  }
}

size_t
SurfaceWaveDispersion::get_number_of_modes_love()
{
  return this->nmodes_love;
}

size_t
SurfaceWaveDispersion::get_number_of_modes_rayleigh()
{
  return this->nmodes_rayleigh;
}

void
SurfaceWaveDispersion::set_number_of_grids(int const new_ngrid)
{
  this->ngrid = new_ngrid;
}

void
SurfaceWaveDispersion::set_basic_parameters(int const new_nf,
                                            double const new_fmax,
                                            double const new_df)
{
  this->nf = new_nf;
  this->fmax = new_fmax;
  this->df = new_df;

  this->generate_frequency_vector();
}

void
SurfaceWaveDispersion::dump_basic_parameters()
{
  std::cout << std::setprecision(8);
  std::cout << "---" << std::endl;
  std::cout << "- Basic params" << std::endl;
  std::cout << "nf = " << this->nf << " points" << std::endl;
  std::cout << "fmax = " << this->fmax << " Hz" << std::endl;
  std::cout << "df = " << this->df << " Hz" << std::endl;
  std::cout << "---" << std::endl;
  std::cout << std::endl;
}

void
SurfaceWaveDispersion::dump_ground_structure()
{
  std::cout << std::endl << "- Thickness" << std::endl;
  for (auto v : thickness) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;

  std::cout << std::endl << "- Interface depth" << std::endl;
  for (auto v : depth) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;

  std::cout << std::endl << "- Vp" << std::endl;
  for (auto v : vp) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;

  std::cout << std::endl << "- Vs" << std::endl;
  for (auto v : vs) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;

  std::cout << std::endl << "- Density" << std::endl;
  for (auto v : rho) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;
  std::cout << std::endl;
}

void
SurfaceWaveDispersion::detect_LVLs()
{
  this->LVL_indices.clear();
  for (size_t i = 1; i < this->vs.size() - 1; ++i) {
    if (this->vs[i] < this->vs[i - 1] && this->vs[i] < this->vs[i + 1]) {
      this->LVL_indices.push_back(i);
    }
  }
}

void
SurfaceWaveDispersion::dump_LVL_indices()
{
  if (not this->LVL_indices.empty()) {
    std::cout << "Indices for LVLs: " << std::endl;
    for (size_t i = 0; i < this->LVL_indices.size(); ++i) {
      std::cout << this->LVL_indices[i] << std::endl;
    }
  } else {
    std::cout << "There is no low velocity layer (LVL)." << std::endl;
  }
  std::cout << std::endl;
}

void
SurfaceWaveDispersion::dump_number_of_modes_love()
{
  std::cout << "---" << std::endl;
  std::cout << "- Number of modes (Love): " << this->nmodes_love << std::endl;
  std::cout << "---" << std::endl;
  std::cout << std::endl;
}

void
SurfaceWaveDispersion::dump_number_of_modes_rayleigh()
{
  std::cout << "---" << std::endl;
  std::cout << "- Number of modes (Rayl): " << this->nmodes_rayleigh
            << std::endl;
  std::cout << "---" << std::endl;
  std::cout << std::endl;
}

void
SurfaceWaveDispersion::inspect_nmodes_love()
{
  for (auto const& row : this->phvs_love) {
    if (row.size() > this->nmodes_love) {
      this->nmodes_love = row.size();
    }
  }
}

void
SurfaceWaveDispersion::inspect_nmodes_rayleigh()
{
  for (auto const& row : this->phvs_rayleigh) {
    if (row.size() > this->nmodes_rayleigh) {
      this->nmodes_rayleigh = row.size();
    }
  }
}

void
SurfaceWaveDispersion::write_love_results(std::string const file_path)
{
  std::ofstream ouf{ file_path };
  ouf << std::setprecision(18) << std::scientific;
  for (size_t i = 0; i < this->phvs_love.size(); ++i) {
    ouf << this->freq[i] << ", ";
    for (size_t j = 0; j < this->nmodes_love - 1; ++j) {
      if (j < this->phvs_love[i].size()) {
        ouf << this->phvs_love[i][j] << ", ";
      } else {
        ouf << std::nan("") << ", ";
      }
    }
    if (this->nmodes_love - 1 < this->phvs_love[i].size()) {
      ouf << this->phvs_love[i][nmodes_love - 1] << std::endl;
    } else {
      ouf << std::nan("") << std::endl;
    }
  }
  ouf.close();
}

void
SurfaceWaveDispersion::write_rayleigh_results(std::string const file_path)
{
  this->write_table(file_path, this->phvs_rayleigh, this->nmodes_rayleigh);
}

void
SurfaceWaveDispersion::write_rayleigh_power_fraction(
  std::string const file_path)
{
  this->write_table(file_path, this->psds_rayleigh, this->nmodes_rayleigh);
}

void
SurfaceWaveDispersion::write_table(
  std::string const file_path,
  std::vector<std::vector<double>> const& table,
  size_t const nmodes)
{
  std::ofstream ouf{ file_path };
  ouf << std::setprecision(18) << std::scientific;
  for (size_t i = 0; i < table.size(); ++i) {
    ouf << this->freq[i] << ", ";
    for (size_t j = 0; j < nmodes - 1; ++j) {
      if (j < table[i].size()) {
        ouf << table[i][j] << ", ";
      } else {
        ouf << std::nan("") << ", ";
      }
    }
    if (nmodes - 1 < table[i].size()) {
      ouf << table[i][nmodes - 1] << std::endl;
    } else {
      ouf << std::nan("") << std::endl;
    }
  }
  ouf.close();
}

void
SurfaceWaveDispersion::read_profile(std::string const path,
                                    char const delimiter)
{
  // delimiter is given as a default argument. -> See header file.
  this->thickness.clear();
  this->vp.clear();
  this->vs.clear();
  this->rho.clear();

  std::ifstream inf{ path };
  std::string line, val;
  while (std::getline(inf, line)) {
    std::stringstream iss(line);
    std::getline(iss, val, delimiter);
    if (val[0] != '#') {
      if (val == "infty") {
        thickness.push_back(std::numeric_limits<double>::infinity());
      } else {
        thickness.push_back(std::stod(val));
      }
      std::getline(iss, val, delimiter);
      this->vp.push_back(std::stod(val));

      std::getline(iss, val, delimiter);
      this->vs.push_back(std::stod(val));

      std::getline(iss, val, delimiter);
      this->rho.push_back(std::stod(val));
    }
  }

  this->thickness2depth();
  this->calc_mu();
  this->detect_LVLs();
}

void
SurfaceWaveDispersion::calc_mu()
{
  this->mu.clear();
  this->mu.resize(this->vs.size());
  for (size_t i = 0; i < this->vs.size(); ++i) {
    this->mu[i] = std::pow(this->vs[i], 2) * this->rho[i];
  }
}

void
SurfaceWaveDispersion::thickness2depth()
{
  this->depth.clear();
  this->depth.resize(this->thickness.size());
  this->depth[0] = 0.0;
  for (size_t i = 0; i < this->thickness.size() - 1; ++i) {
    this->depth[i + 1] = this->depth[i] + this->thickness[i];
  }
  this->depth.push_back(std::numeric_limits<double>::infinity());
}

double
SurfaceWaveDispersion::vector_rms(std::vector<double> const& vec)
{
  double ans = 0.0;
  for (auto v : vec) {
    ans += std::pow(v, 2);
  }
  return std::sqrt(ans / (double)vec.size());
}

double
SurfaceWaveDispersion::vector_rms(std::vector<double> const& vec,
                                  size_t ind_s,
                                  size_t ind_e)
{
  double ans = 0.0;
  for (size_t i = ind_s; i <= ind_e; ++i) {
    ans += std::pow(vec[i], 2);
  }
  return std::sqrt(ans / (double)(ind_e - ind_s + 1));
}

// bisect Rev.
double
SurfaceWaveDispersion::bisection_sf(double const x_min,
                                    double const x_max,
                                    double const omega,
                                    int const key,
                                    std::string target)

{
  std::function<std::vector<double>(double, double)> target_function;
  if (target == "rayleigh") {
    target_function = [this](double phv, double omega) {
      return this->rayleigh_secular_function(phv, omega);
    };
  } else {
    perror("Error -> Bisection No. 0");
    std::cerr << "Target should be either \"love\" or \"rayleigh\"."
              << std::endl;
    exit(errno);
  }

  // key for rms calculation
  int key_1, key_2;
  if (key % 2 == 0) {
    key_1 = key;
    key_2 = key + 1;
  } else {
    key_1 = key - 1;
    key_2 = key;
  }

  double y_l = target_function(x_min, omega)[key];
  double y_r = target_function(x_max, omega)[key];
  double y_tot_l = vector_rms(target_function(x_min, omega), key_1, key_2);
  double y_tot_r = vector_rms(target_function(x_max, omega), key_1, key_2);
  double x_l = x_min;
  double x_r = x_max;

  if (y_l * y_r > 0.0) {
    perror("Error -> Bisection No. 1");
    exit(errno);
  }

  double const ytol =
    std::abs(y_tot_l - y_tot_r) / (x_r - x_l) * this->bisect_acc * 0.1;

  for (int i = 0; i < this->bisect_max_itr; ++i) {
    double x_c = (x_l + x_r) / 2.0;

    if (x_r - x_l < this->bisect_acc) {
      if (vector_rms(target_function(x_c, omega), key_1, key_2) < ytol) {
        return x_c;
      }
    }

    if (i == this->bisect_max_itr - 1) {
      if (not(vector_rms(target_function(x_c, omega), key_1, key_2) < ytol)) {
        return std::nan("");
      }
      perror("Error -> Bisection No. 2");
      exit(errno);
    }

    double y_c = target_function(x_c, omega)[key];
    y_l = target_function(x_l, omega)[key];
    if (y_c * y_l <= 0.0) {
      x_r = x_c;
    } else {
      x_l = x_c;
    }
  }

  perror("Error -> Bisection No. 3");
  exit(errno);
}

std::complex<double>
SurfaceWaveDispersion::love_secular_function(double const phv,
                                             double const omega)
{
  std::complex<double> k = omega / phv;
  std::complex<double> imunit{ 0.0, 1.0 };

  std::complex<double> ans = std::nan("");

  std::complex<double> R_du_ip1{ 0.0, 0.0 }; // \hat{R_{du}}^{(N + 1)} = 0
  for (int i = (int)this->thickness.size() - 2; i >= 0; --i) {
    std::complex<double> nu_i = this->nu(k, omega, i);
    std::complex<double> nu_ip1 = this->nu(k, omega, i + 1);

    std::complex<double> mu_i = this->mu[i];
    std::complex<double> mu_ip1 = this->mu[i + 1];

    std::complex<double> Lambda_u_ip1_zi =
      std::exp(-nu_ip1 * this->thickness[i + 1]);

    std::complex<double> Lambda_d_i_zi = std::exp(-nu_i * this->thickness[i]);

    std::complex<double> a = 0.5 * (1.0 + mu_ip1 * nu_ip1 / (mu_i * nu_i));
    std::complex<double> b = 0.5 * (1.0 - mu_ip1 * nu_ip1 / (mu_i * nu_i));

    std::complex<double> T_d_i =
      Lambda_d_i_zi / (a + b * Lambda_u_ip1_zi * R_du_ip1);
    std::complex<double> R_du_i = (b + a * Lambda_u_ip1_zi * R_du_ip1) * T_d_i;
    if (i == 0) {
      std::complex<double> E21 = -mu_i * nu_i;
      std::complex<double> E22 = mu_i * nu_i;
      std::complex<double> Lambda_u_1_z0 = std::exp(-nu_i * this->thickness[i]);
      ans = (E21 + E22 * Lambda_u_1_z0 * R_du_i);
      break;
    }
    R_du_ip1 = R_du_i;
  }
  return ans;
}

std::complex<double>
SurfaceWaveDispersion::nu(std::complex<double> const k,
                          double const omega,
                          int const ind)
{
  std::complex<double> ans =
    std::sqrt(k * k - std::pow(omega / this->vs[ind], 2));
  if (ans.real() < 0) {
    ans *= -1.0;
  }
  return ans;
}

std::complex<double>
SurfaceWaveDispersion::gamma(std::complex<double> const k,
                             double const omega,
                             int const ind)
{
  std::complex<double> ans =
    std::sqrt(k * k - std::pow(omega / this->vp[ind], 2));
  if (ans.real() < 0) {
    ans *= -1.0;
  }
  return ans;
}

std::complex<double>
SurfaceWaveDispersion::chi1(std::complex<double> const k,
                            double const omega,
                            int const ind)
{

  std::complex<double> ans = 2.0 * k - std::pow(omega / this->vs[ind], 2) / k;
  return ans;
}

Eigen::Matrix<std::complex<double>, 2, 2>
SurfaceWaveDispersion::Lambda_d(int const ind,
                                double const z,
                                std::complex<double> const gamma_ind,
                                std::complex<double> const nu_ind)
{
  Eigen::Matrix<std::complex<double>, 2, 2> ans =
    Eigen::Matrix<std::complex<double>, 2, 2>::Zero();
  ans(0, 0) = std::exp(-gamma_ind * (z - this->depth[ind]));
  ans(1, 1) = std::exp(-nu_ind * (z - this->depth[ind]));
  // z_{i - 1} = depth[i]
  return ans;
}

Eigen::Matrix<std::complex<double>, 2, 2>
SurfaceWaveDispersion::Lambda_u(int const ind,
                                double const z,
                                std::complex<double> const gamma_ind,
                                std::complex<double> const nu_ind)
{
  Eigen::Matrix<std::complex<double>, 2, 2> ans =
    Eigen::Matrix<std::complex<double>, 2, 2>::Zero();
  ans(0, 0) = std::exp(-gamma_ind * (this->depth[ind + 1] - z));
  ans(1, 1) = std::exp(-nu_ind * (this->depth[ind + 1] - z));
  // z_{i} = depth[i + 1]
  return ans;
}

bool
SurfaceWaveDispersion::is_LVL(int const ind)
{
  bool ans = false;
  for (int i = 0; i < (int)this->LVL_indices.size(); ++i) {
    if (this->LVL_indices[i] == ind) {
      return true;
    }
  }
  return ans;
}

void
SurfaceWaveDispersion::rayleigh_downward_recursion(
  std::complex<double> const k,
  double const omega,
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>>& vec_E_i,
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>>& vec_E_inv_i,
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>>& vec_Lambda_kernel_i,
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>>& vec_R_du_i)
{
  std::complex<double> const one{ 1.0, 0.0 };

  int const n_layers = this->thickness.size(); // N + 1

  vec_E_i.resize(n_layers);
  vec_E_inv_i.resize(n_layers);
  vec_Lambda_kernel_i.resize(n_layers);
  vec_R_du_i.resize(n_layers); // R_{du}^{(N + 1)} to R_{du}^{(1)}

  // for ind = 1 to N+1 (i.e., for each layer)
  for (int i = 0; i < n_layers; ++i) {
    std::complex<double> const mu_i = this->mu[i];
    std::complex<double> const nu_i = this->nu(k, omega, i);
    std::complex<double> const gamma_i = this->gamma(k, omega, i);
    std::complex<double> const chi1_i = this->chi1(k, omega, i);

    std::complex<double> const nu1_i = nu_i / k;
    std::complex<double> const gamma1_i = gamma_i / k;

    vec_Lambda_kernel_i[i] = Eigen::Matrix<std::complex<double>, 2, 2>::Zero();
    vec_Lambda_kernel_i[i](0, 0) = std::exp(-gamma_i * this->thickness[i]);
    vec_Lambda_kernel_i[i](1, 1) = std::exp(-nu_i * this->thickness[i]);

    // clang-format off
    vec_E_inv_i[i] <<
      -2.0*k*mu_i,       mu_i*chi1_i/gamma1_i,  one/gamma1_i,  -one,
      -mu_i*chi1_i/nu1_i, 2.0*k*mu_i,           one,           -one/nu1_i,
      -2.0*k*mu_i,       -mu_i*chi1_i/gamma1_i, -one/gamma1_i, -one,
      -mu_i*chi1_i/nu1_i, -2.0*k*mu_i,          -one,          -one/nu1_i;
    vec_E_inv_i[i] /= (mu_i * (4.0 * k - 2.0 * chi1_i));
    // clang-format on

    // clang-format off
    vec_E_i[i] <<
      -one,             nu1_i,          -one,              nu1_i,
      -gamma1_i,        one,            gamma1_i,          -one,
      2.0*mu_i*gamma_i, -chi1_i*mu_i,   -2.0*mu_i*gamma_i, chi1_i*mu_i,
      chi1_i*mu_i,      -2.0*mu_i*nu_i, chi1_i*mu_i,       -2.0*mu_i*nu_i;
    // clang-format on
  }

  // \hat{R_{du}}^{(N + 1)} = 0 >>> Initial point
  vec_R_du_i[n_layers - 1] = Eigen::Matrix<std::complex<double>, 2, 2>::Zero();

  // Propagate matrices (downgoing)
  // for ind = N to 1
  // (ind = i + 1) -> i = N - 1 to 0
  for (int i = n_layers - 2; i >= 0; --i) {
    // \Lambda_{u}^{(i + 1)}(z^{i})
    Eigen::Matrix<std::complex<double>, 2, 2> Lambda_u_ip1_zi =
      vec_Lambda_kernel_i[i + 1];

    // \Lambda_{d}^{(i)}(z^{i})
    Eigen::Matrix<std::complex<double>, 2, 2> Lambda_d_i_zi =
      vec_Lambda_kernel_i[i];

    // (E^{(i)})^{-1}
    Eigen::Matrix<std::complex<double>, 4, 4> E_inv_i = vec_E_inv_i[i];

    // E^{(i + 1)}
    Eigen::Matrix<std::complex<double>, 4, 4> E_ip1 = vec_E_i[i + 1];

    // W^{(i)}
    Eigen::Matrix<std::complex<double>, 4, 4> W = E_inv_i * E_ip1;
    Eigen::Matrix<std::complex<double>, 2, 2> W11 = W.block<2, 2>(0, 0);
    Eigen::Matrix<std::complex<double>, 2, 2> W12 = W.block<2, 2>(0, 2);
    Eigen::Matrix<std::complex<double>, 2, 2> W21 = W.block<2, 2>(2, 0);
    Eigen::Matrix<std::complex<double>, 2, 2> W22 = W.block<2, 2>(2, 2);

    // \hat{T}_d^{i}
    Eigen::Matrix<std::complex<double>, 2, 2> T_d_i =
      (W11 + W12 * Lambda_u_ip1_zi * vec_R_du_i[i + 1]).inverse() *
      Lambda_d_i_zi;

    // \hat{R}_{du}^{i}
    vec_R_du_i[i] = (W21 + W22 * Lambda_u_ip1_zi * vec_R_du_i[i + 1]) * T_d_i;
  }
}

Eigen::Matrix<std::complex<double>, 2, 2>
SurfaceWaveDispersion::rayleigh_surface_response(std::complex<double> const k,
                                                 double const omega)
{
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>> vec_E_i;
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>> vec_E_inv_i;
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>> vec_Lambda_kernel_i;
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>> vec_R_du_i;

  this->rayleigh_downward_recursion(
    k, omega, vec_E_i, vec_E_inv_i, vec_Lambda_kernel_i, vec_R_du_i);

  // The wavefield in the top layer is expressed by its downgoing
  // coefficients C_d only (C_u = \hat{R}_{du}^{(1)} C_d). With
  // L = \Lambda_u^{(1)}(z^{(0)}) \hat{R}_{du}^{(1)}, at the free surface
  //   displacement = (E11 + E12 L) C_d,
  //   traction = (E21 + E22 L) C_d.
  // The determinant of the latter matrix is the secular function F_0.
  Eigen::Matrix<std::complex<double>, 4, 4> const E_1 = vec_E_i[0];
  Eigen::Matrix<std::complex<double>, 2, 2> const E11 = E_1.block<2, 2>(0, 0);
  Eigen::Matrix<std::complex<double>, 2, 2> const E12 = E_1.block<2, 2>(0, 2);
  Eigen::Matrix<std::complex<double>, 2, 2> const E21 = E_1.block<2, 2>(2, 0);
  Eigen::Matrix<std::complex<double>, 2, 2> const E22 = E_1.block<2, 2>(2, 2);

  Eigen::Matrix<std::complex<double>, 2, 2> const Lambda_R =
    vec_Lambda_kernel_i[0] * vec_R_du_i[0];

  Eigen::Matrix<std::complex<double>, 2, 2> const D = E11 + E12 * Lambda_R;
  Eigen::Matrix<std::complex<double>, 2, 2> const S = E21 + E22 * Lambda_R;

  return D * S.inverse();
}

std::complex<double>
SurfaceWaveDispersion::rayleigh_medium_response(double const phv,
                                                double const omega,
                                                double const radius)
{
  // Res_{k = k_n} Y_zz(k) = 1 / (2 pi i) * \oint Y_zz(k) dk on the circle
  // k = k_n + radius * exp(i theta), evaluated by the trapezoidal rule. Y is
  // analytic around the circle (it is even in nu and gamma of every finite
  // layer), which makes the rule spectrally accurate and, unlike a finite
  // difference of the secular function, insensitive to the step size and to
  // a small error of the root.
  std::complex<double> const imunit{ 0.0, 1.0 };
  double const k_n = omega / phv;

  std::complex<double> residue{ 0.0, 0.0 };
  for (int m = 0; m < this->contour_npts; ++m) {
    double const theta =
      2.0 * std::numbers::pi * ((double)m + 0.5) / (double)this->contour_npts;
    std::complex<double> const dk = radius * std::exp(imunit * theta);
    residue += dk * this->rayleigh_surface_response(k_n + dk, omega)(1, 1);
  }
  residue /= (double)this->contour_npts;

  return k_n * residue;
}

std::vector<double>
SurfaceWaveDispersion::rayleigh_secular_function(double const phv,
                                                 double const omega)
{
  std::complex<double> const k{ omega / phv, 0.0 };
  std::complex<double> const imunit{ 0.0, 1.0 };
  std::complex<double> const zero{ 0.0, 0.0 };
  Eigen::Matrix<std::complex<double>, 2, 2> const Imat_2x2 =
    Eigen::Matrix<std::complex<double>, 2, 2>::Identity();

  std::vector<double> ans;

  int const n_layers = this->thickness.size(); // N + 1

  // Generalized R/T coefficients
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>> vec_R_du_i;
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>> vec_R_ud_i;

  vec_R_ud_i.resize(n_layers - 1); // R_{ud}^{(0)} to R_{ud}^{(N - 1)}
  // Note that
  // vec_R_du_i[i] = R_{du}^{(i + 1)} for i = 0, 1, ..., N
  // but
  // vec_R_ud_i[i] = R_{ud}^{(i)} for i = 0, 1, ..., N.

  // Matrices for each layer.
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>> vec_E_i;
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>> vec_E_inv_i;
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>> vec_Lambda_kernel_i;

  // Propagate matrices (downgoing)
  this->rayleigh_downward_recursion(
    k, omega, vec_E_i, vec_E_inv_i, vec_Lambda_kernel_i, vec_R_du_i);

  if (n_layers >= 2) {
    Eigen::Matrix<std::complex<double>, 4, 4> E_i = vec_E_i[0];
    Eigen::Matrix<std::complex<double>, 2, 2> E21 = E_i.block<2, 2>(2, 0);
    Eigen::Matrix<std::complex<double>, 2, 2> E22 = E_i.block<2, 2>(2, 2);

    Eigen::Matrix<std::complex<double>, 2, 2> Lambda_u_1_z0 =
      vec_Lambda_kernel_i[0];

    std::complex<double> r_scf0 =
      (E21 + E22 * Lambda_u_1_z0 * vec_R_du_i[0]).determinant();
    ans.push_back(r_scf0.real());
    ans.push_back(r_scf0.imag());

    vec_R_ud_i[0] = -E21.inverse() * E22 * Lambda_u_1_z0;

    // std::complex<double> r_scf1 =
    //   (Imat_2x2 - vec_R_ud_i[0] * vec_R_du_i[0]).determinant();
    // ans.push_back(r_scf1.real());
    // ans.push_back(r_scf1.imag());
  }

  if (this->LVL_indices.empty()) {
    return ans;
  }

  // Propagate matrices (upgoing)
  // for ind = 1 to j_{max}
  // (ind = i) -> i = 1 to j_{max}

  int const max_LVL_ind = this->LVL_indices[this->LVL_indices.size() - 1];
  for (int i = 1; i < max_LVL_ind + 1; ++i) {
    // \Lambda_{u}^{(i + 1)}(z^{i})
    Eigen::Matrix<std::complex<double>, 2, 2> Lambda_u_ip1_zi =
      vec_Lambda_kernel_i[i];

    // \Lambda_{d}^{(i)}(z^{i})
    Eigen::Matrix<std::complex<double>, 2, 2> Lambda_d_i_zi =
      vec_Lambda_kernel_i[i - 1];

    // (E^{(i)})^{-1}
    Eigen::Matrix<std::complex<double>, 4, 4> E_inv_ip1 = vec_E_inv_i[i];

    // E^{(i + 1)}
    Eigen::Matrix<std::complex<double>, 4, 4> E_i = vec_E_i[i - 1];

    // W^{(i)}
    Eigen::Matrix<std::complex<double>, 4, 4> O = E_inv_ip1 * E_i;
    Eigen::Matrix<std::complex<double>, 2, 2> O11 = O.block<2, 2>(0, 0);
    Eigen::Matrix<std::complex<double>, 2, 2> O12 = O.block<2, 2>(0, 2);
    Eigen::Matrix<std::complex<double>, 2, 2> O21 = O.block<2, 2>(2, 0);
    Eigen::Matrix<std::complex<double>, 2, 2> O22 = O.block<2, 2>(2, 2);

    // \hat{T}_d^{i}
    Eigen::Matrix<std::complex<double>, 2, 2> T_u_i =
      (O21 * Lambda_d_i_zi * vec_R_ud_i[i - 1] + O22).inverse() *
      Lambda_u_ip1_zi;

    // \hat{R}_{du}^{i}
    vec_R_ud_i[i] = (O11 * Lambda_d_i_zi * vec_R_ud_i[i - 1] + O12) * T_u_i;

    // if (i - 1) in ind list of LVL
    // IMPL: change here
    if (this->is_LVL(i)) {
      std::complex<double> r_scfj =
        (Imat_2x2 - vec_R_ud_i[i] * vec_R_du_i[i]).determinant();
      ans.push_back(r_scfj.real());
      ans.push_back(r_scfj.imag());
    }
  }
  return ans;
}

void
SurfaceWaveDispersion::calculate_love_power_fraction()
{
  std::cerr << "calculate_love_power_fraction() has not been implemented yet!"
            << std::endl;
}

void
SurfaceWaveDispersion::calculate_rayleigh_power_fraction()
{
  // Power of each mode relative to that of the 0th mode, following the
  // Thomson-Haskell implementation: (A_j^2 c_j) / (A_0^2 c_0), where A is the
  // medium response (Harkrider, 1964) for a vertical point force and the
  // vertical displacement at the free surface. Here A is evaluated as the
  // residue of the surface response; see rayleigh_medium_response().
  this->psds_rayleigh.clear();

  size_t n_unresolved = 0;
  size_t n_invalid_0th = 0;

  for (size_t i = 0; i < this->phvs_rayleigh.size(); ++i) {
    std::vector<double> const& phvs_at_freq = this->phvs_rayleigh[i];
    std::vector<double> psds_at_freq(phvs_at_freq.size(), std::nan(""));

    // freq[0] = 0: the medium response is not defined.
    if (i == 0) {
      this->psds_rayleigh.push_back(psds_at_freq);
      continue;
    }

    double const omega = this->freq[i] * 2.0 * std::numbers::pi;

    double power_0th = std::nan("");
    for (size_t j = 0; j < phvs_at_freq.size(); ++j) {
      double const phv = phvs_at_freq[j];
      double const k_n = omega / phv;

      // The circle must not enclose the other roots.
      double radius = this->contour_rel_radius * k_n;
      for (size_t l = 0; l < phvs_at_freq.size(); ++l) {
        if (l != j) {
          radius = std::min(radius,
                            this->contour_gap_fraction *
                              std::abs(omega / phvs_at_freq[l] - k_n));
        }
      }

      // The circle must enclose the true root, whose position is known only
      // within bisect_acc / 2 in phase velocity. Give up if two roots are too
      // close to each other to guarantee that.
      double const k_uncertainty = k_n * 0.5 * this->bisect_acc / phv;
      if (not(radius > 4.0 * k_uncertainty)) {
        ++n_unresolved;
        continue;
      }

      std::complex<double> const medium_response =
        this->rayleigh_medium_response(phv, omega, radius);
      double const power = std::norm(medium_response) * phv;

      if (j == 0) {
        power_0th = power;
      }
      psds_at_freq[j] = power / power_0th;
    }

    if (not phvs_at_freq.empty() &&
        not(std::isfinite(power_0th) && power_0th > 0.0)) {
      ++n_invalid_0th;
    }
    this->psds_rayleigh.push_back(psds_at_freq);
  }

  if (n_unresolved > 0) {
    std::cerr << "[Rayl power fraction] " << n_unresolved
              << " mode(s) were too close to a neighboring root to be "
                 "resolved (written as nan)."
              << std::endl;
  }
  if (n_invalid_0th > 0) {
    std::cerr << "[Rayl power fraction] The power of the 0th mode was zero or "
                 "not available at "
              << n_invalid_0th
              << " frequency point(s); the ratios there are not meaningful."
              << std::endl;
  }
}

double
SurfaceWaveDispersion::predicted_num_of_modes(double phv, double omega)
{
  double ans = 0.0;
  for (size_t i = 0; i < this->thickness.size() - 2; ++i) {
    std::complex<double> gam_p = std::sqrt(
      std::complex<double>{ std::pow(phv / this->vp[i], 2) - 1.0, 0 });
    // if (gam_p.real() < 0) {
    //   gam_p *= -1.0;
    // }

    std::complex<double> gam_s = std::sqrt(
      std::complex<double>{ std::pow(phv / this->vs[i], 2) - 1.0, 0 });
    // if (gam_s.real() < 0) {
    //   gam_s *= -1.0;
    // }
    // std::cerr << "gam_p, gam_s = " << gam_p << ", " << gam_s << std::endl;

    ans += ((gam_s + gam_p) * this->thickness[i]).real();
  }
  // ans = omega / (2 * pi) * 2 * (sum (gamma + nu)*h) / phv
  ans *= omega / (std::numbers::pi * phv);

  return ans;
}

bool
SurfaceWaveDispersion::equals_vs_of_a_layer(double phv)
{
  bool ans = false;

  for (size_t i = 0; i < this->vs.size(); ++i) {
    // if (abs(phv - this->vs[i]) < this->bisect_acc) {
    if (abs(phv - this->vs[i]) < 1e-3) {
      ans = true;
    }
  }

  return ans;
}

bool
SurfaceWaveDispersion::equals_vp_of_a_layer(double phv)
{
  bool ans = false;

  for (size_t i = 0; i < this->vp.size(); ++i) {
    // if (abs(phv - this->vp[i]) < this->bisect_acc) {
    if (abs(phv - this->vp[i]) < 1e-3) {
      ans = true;
    }
  }

  return ans;
}

void
SurfaceWaveDispersion::calculate_rayleigh_phase_velocities()
{
  this->phvs_rayleigh.clear();
  double const vs_half_space = this->vs[this->vs.size() - 1];

  phvs_rayleigh.push_back(std::vector<double>{ this->vs[this->vs.size() - 1] });
  for (size_t i = 1; i < this->freq.size(); ++i) {
    std::vector<double> phvs_rayleigh_at_freq;

    double omega = freq[i] * 2.0 * std::numbers::pi;

    // Predicte number of modes using Fan et al. (2007), published by Chinese
    // Journal of Geophysics. Note that it depends on frequency.
    int const num_of_grid_pts =
      (this->predicted_num_of_modes(vs_half_space, omega) + 1) * this->ngrid;

    double min_phv_bf = 0.0;
    if (i > 1) {
      min_phv_bf = *std::min_element(this->phvs_rayleigh[i - 1].begin(),
                                     this->phvs_rayleigh[i - 1].end());
    }
    double const min_vs = *std::min_element(this->vs.begin(), this->vs.end());
    double const dc = (vs_half_space - 0.8 * min_vs) / (double)num_of_grid_pts;

    // Prepare grid points for grid search
    std::vector<double> grid_pts;
    grid_pts.resize(num_of_grid_pts);
    for (int j = 0; j < num_of_grid_pts; ++j) {
      grid_pts[j] =
        0.8 * min_vs + (double)(j + 1) * dc + this->fluctuation_for_root_search;
    }
    // if (i != 0) {
    if (i > 1) {
      for (size_t j = 0; j < phvs_rayleigh[i - 1].size(); ++j) {
        double center = phvs_rayleigh[i - 1][j];
        // double center = 2.0 * phvs_rayleigh[i - 1][j] - phvs_rayleigh[i -
        // 2][j];

        //  Following an idea presented by Fan et al. (2007).
        // std::vector<double> additional_grid_pts{ -5.0, -0.5, 0, 0.5, 5.0 };
        std::vector<double> additional_grid_pts{
          -5.0,    -0.5,   -0.05, -0.005, -0.0005, -0.00005, 0,
          0.00005, 0.0005, 0.005, 0.05,   0.5,     5.0
        };
        for (size_t k = 0; k < additional_grid_pts.size(); ++k) {
          additional_grid_pts[k] += center;
          grid_pts.push_back(additional_grid_pts[k]);
        }
      }
    }

    std::sort(grid_pts.begin(), grid_pts.end());

    // Grid search
    double phv_bf = grid_pts[0];
    std::vector<double> bf = rayleigh_secular_function(phv_bf, omega);
    for (size_t j = 1; j < grid_pts.size(); ++j) {
      double phv_af = grid_pts[j];
      std::vector<double> af = rayleigh_secular_function(phv_af, omega);

      // Check kth secular function
      for (size_t k = 0; k < af.size(); ++k) {
        if (bf[k] * af[k] < 0 && std::abs(bf[k] * af[k]) > 1e-8) {
          double candidate = bisection_sf(phv_bf, phv_af, omega, k, "rayleigh");
          if (not std::isnan(candidate) && not std::isinf(candidate) &&
              not this->equals_vs_of_a_layer(candidate) &&
              not this->equals_vp_of_a_layer(candidate)) {
            phvs_rayleigh_at_freq.push_back(candidate);
            break;
          }
        }
      }
      phv_bf = phv_af;
      bf = std::move(af);
    }
    // ouf.close();
    this->phvs_rayleigh.push_back(phvs_rayleigh_at_freq);
  }
  if (phvs_rayleigh[1].size() > 0 and phvs_rayleigh[2].size() > 0) {
    this->phvs_rayleigh[0][0] =
      this->phvs_rayleigh[1][0] * 2 - this->phvs_rayleigh[2][0];
  } else {
    std::cerr << "Something went wrong in low frequency region." << std::endl;
    std::cerr << "Increasing ngrid might resolve the problem." << std::endl;
    this->phvs_rayleigh[0][0] = this->phvs_rayleigh[1][0];
  }
  inspect_nmodes_rayleigh();
}
