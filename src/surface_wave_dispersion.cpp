// SPDX-FileCopyrightText: Copyright 2024, 2026 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "surface_wave_dispersion.hpp"

SurfaceWaveDispersion::SurfaceWaveDispersion(int const new_nf,
                                             double const new_fmin,
                                             double const new_fmax)
{
  this->set_basic_parameters(new_nf, new_fmin, new_fmax);
}

void
SurfaceWaveDispersion::generate_frequency_vector()
{
  this->freq.clear();
  this->freq.resize(this->nf + 1);
  for (int i = 0; i < this->nf + 1; ++i) {
    this->freq[i] =
      this->fmin + (this->fmax - this->fmin) / (double)this->nf * (double)i;
  }
}

size_t
SurfaceWaveDispersion::get_number_of_modes_rayleigh()
{
  return this->nmodes_rayleigh;
}

bool
SurfaceWaveDispersion::has_attenuation()
{
  return this->attenuation_given;
}

void
SurfaceWaveDispersion::set_ignore_attenuation(
  bool const new_attenuation_ignored)
{
  this->attenuation_ignored = new_attenuation_ignored;
}

void
SurfaceWaveDispersion::set_number_of_grids(int const new_ngrid)
{
  this->ngrid = new_ngrid;
}

void
SurfaceWaveDispersion::set_basic_parameters(int const new_nf,
                                            double const new_fmin,
                                            double const new_fmax)
{
  this->nf = new_nf;
  this->fmin = new_fmin;
  this->fmax = new_fmax;
  this->df = (new_fmax - new_fmin) / (double)new_nf;

  this->generate_frequency_vector();
}

void
SurfaceWaveDispersion::set_reference_frequency(double const new_f_ref)
{
  this->reference_frequency_given = true;
  this->f_ref = new_f_ref;
}

void
SurfaceWaveDispersion::dump_basic_parameters()
{
  std::cout << std::setprecision(8);
  std::cout << "---" << std::endl;
  std::cout << "- Basic params" << std::endl;
  std::cout << "nf = " << this->nf << " points" << std::endl;
  std::cout << "fmin = " << this->fmin << " Hz" << std::endl;
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
  for (auto v : vp_ref) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;

  std::cout << std::endl << "- Vs" << std::endl;
  for (auto v : vs_ref) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;

  std::cout << std::endl << "- Density" << std::endl;
  for (auto v : rho) {
    std::cout << v << std::endl;
  }
  std::cout << "---" << std::endl;

  if (this->attenuation_given) {
    std::cout << std::endl << "- Qp" << std::endl;
    for (auto v : qp) {
      std::cout << v << std::endl;
    }
    std::cout << "---" << std::endl;

    std::cout << std::endl << "- Qs" << std::endl;
    for (auto v : qs) {
      std::cout << v << std::endl;
    }
    std::cout << "---" << std::endl;

    if (this->reference_frequency_given) {
      std::cout
        << std::endl
        << "- Velocity model: Kjartansson (1979), reference frequency = "
        << this->f_ref << " Hz" << std::endl;
    } else {
      std::cout << std::endl
                << "- Velocity model: frequency independent" << std::endl;
    }
    std::cout << "---" << std::endl;
  }
  std::cout << std::endl;
}

std::vector<int>
SurfaceWaveDispersion::detect_LVLs(std::vector<std::complex<double>> const& vs)
{
  std::vector<int> ans;
  for (size_t i = 1; i + 1 < vs.size(); ++i) {
    if (vs[i].real() < vs[i - 1].real() && vs[i].real() < vs[i + 1].real()) {
      ans.push_back(i);
    }
  }
  return ans;
}

void
SurfaceWaveDispersion::dump_LVL_indices()
{
  std::vector<int> const& LVL_indices = this->medium_elastic_ref.LVL_indices;
  if (not LVL_indices.empty()) {
    if (this->reference_frequency_given) {
      std::cout << "Indices for LVLs (at the reference frequency): "
                << std::endl;
    } else {
      std::cout << "Indices for LVLs: " << std::endl;
    }
    for (size_t i = 0; i < LVL_indices.size(); ++i) {
      std::cout << LVL_indices[i] << std::endl;
    }
  } else {
    std::cout << "There is no low velocity layer (LVL)." << std::endl;
  }
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
SurfaceWaveDispersion::inspect_nmodes_rayleigh()
{
  for (auto const& row : this->phvs_rayleigh) {
    if (row.size() > this->nmodes_rayleigh) {
      this->nmodes_rayleigh = row.size();
    }
  }
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
SurfaceWaveDispersion::write_rayleigh_medium_response(
  std::string const file_path)
{
  this->write_table(
    file_path, this->medium_responses_rayleigh, this->nmodes_rayleigh);
}

void
SurfaceWaveDispersion::write_rayleigh_attenuation(std::string const file_path)
{
  this->write_table(file_path, this->attens_rayleigh, this->nmodes_rayleigh);
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
    ouf << this->freq[i];
    for (size_t j = 0; j < nmodes; ++j) {
      // A mode that does not exist is written as nan (without a sign).
      if (j < table[i].size() && not std::isnan(table[i][j])) {
        ouf << ", " << table[i][j];
      } else {
        ouf << ", " << std::nan("");
      }
    }
    ouf << std::endl;
  }
  ouf.close();
}

void
SurfaceWaveDispersion::read_profile(std::string const path,
                                    char const delimiter)
{
  // delimiter is given as a default argument. -> See header file.
  // One line per layer: thickness, Vp, Vs, rho (elastic medium) or
  // thickness, Vp, Vs, rho, Qp, Qs (viscoelastic medium).
  this->thickness.clear();
  this->vp_ref.clear();
  this->vs_ref.clear();
  this->rho.clear();
  this->qp.clear();
  this->qs.clear();

  std::ifstream inf{ path };
  if (not inf) {
    std::cerr << "Error: Cannot open " << path << std::endl;
    exit(EXIT_FAILURE);
  }

  auto const trim = [](std::string const& str) {
    size_t const first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
      return std::string{};
    }
    size_t const last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
  };

  size_t n_columns = 0;
  std::string line, val;
  while (std::getline(inf, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::vector<std::string> fields;
    std::stringstream iss(line);
    while (std::getline(iss, val, delimiter)) {
      fields.push_back(trim(val));
    }

    if (fields.size() != 4 && fields.size() != 6) {
      std::cerr << "Error: A layer must be given by 4 values (thickness, Vp, "
                   "Vs, rho) or 6 values (thickness, Vp, Vs, rho, Qp, Qs): "
                << line << std::endl;
      exit(EXIT_FAILURE);
    }
    if (n_columns == 0) {
      n_columns = fields.size();
    } else if (n_columns != fields.size()) {
      std::cerr << "Error: All the layers must have the same number of values."
                << std::endl;
      exit(EXIT_FAILURE);
    }

    if (fields[0] == "infty") {
      this->thickness.push_back(std::numeric_limits<double>::infinity());
    } else {
      this->thickness.push_back(std::stod(fields[0]));
    }
    this->vp_ref.push_back(std::stod(fields[1]));
    this->vs_ref.push_back(std::stod(fields[2]));
    this->rho.push_back(std::stod(fields[3]));
    if (fields.size() == 6) {
      this->qp.push_back(std::stod(fields[4]));
      this->qs.push_back(std::stod(fields[5]));
    }
  }

  if (this->thickness.empty()) {
    std::cerr << "Error: No layer was found in " << path << std::endl;
    exit(EXIT_FAILURE);
  }

  this->attenuation_given = (n_columns == 6) && not this->attenuation_ignored;

  this->thickness2depth();
  this->build_reference_media();
}

void
SurfaceWaveDispersion::build_reference_media()
{
  // Complex velocity of the constant-Q model of Kjartansson (1979, JGR 84,
  // 4737-4748) at the reference frequency. With gamma = arctan(1 / Q) / pi,
  // the phase velocity c_0 read from the profile and the complex modulus M_0
  // are related by c_0 = sqrt(M_0 / rho) / cos(pi * gamma / 2) (eq. 26), so
  // that the complex velocity sqrt(M_0 / rho) is
  //   c_0 * cos(pi * gamma / 2) * exp(i * pi * gamma / 2)
  //   = c_0 * (1 + s + i / Q) / (2 * s), s = sqrt(1 + Q^{-2}).
  // The time dependence is exp(+i * omega * t), as in the paper (eq. 12).
  // The same complex velocity is used at all the frequencies if the
  // reference frequency is not given.
  auto const complex_factor = [](double const q) {
    if (not(q > 0.0)) {
      return std::complex<double>{ 1.0, 0.0 };
    }
    double const s = std::sqrt(1.0 + 1.0 / (q * q));
    return std::complex<double>{ 1.0 + s, 1.0 / q } / (2.0 * s);
  };

  size_t const n_layers = this->thickness.size();

  Medium& elastic = this->medium_elastic_ref;
  Medium& viscoelastic = this->medium_viscoelastic_ref;
  for (Medium* medium : { &elastic, &viscoelastic }) {
    medium->vp.resize(n_layers);
    medium->vs.resize(n_layers);
    medium->mu.resize(n_layers);
  }
  this->gamma_p.assign(n_layers, 0.0);
  this->gamma_s.assign(n_layers, 0.0);

  for (size_t i = 0; i < n_layers; ++i) {
    elastic.vp[i] = this->vp_ref[i];
    elastic.vs[i] = this->vs_ref[i];
    elastic.mu[i] = std::pow(this->vs_ref[i], 2) * this->rho[i];

    viscoelastic.vp[i] = elastic.vp[i];
    viscoelastic.vs[i] = elastic.vs[i];
    if (this->attenuation_given) {
      viscoelastic.vp[i] = this->vp_ref[i] * complex_factor(this->qp[i]);
      viscoelastic.vs[i] = this->vs_ref[i] * complex_factor(this->qs[i]);
      if (this->qp[i] > 0.0) {
        this->gamma_p[i] = std::atan(1.0 / this->qp[i]) / std::numbers::pi;
      }
      if (this->qs[i] > 0.0) {
        this->gamma_s[i] = std::atan(1.0 / this->qs[i]) / std::numbers::pi;
      }
    }
    viscoelastic.mu[i] = viscoelastic.vs[i] * viscoelastic.vs[i] * this->rho[i];
  }

  elastic.LVL_indices = this->detect_LVLs(elastic.vs);
  viscoelastic.LVL_indices = elastic.LVL_indices;
}

Medium
SurfaceWaveDispersion::medium_at(double const omega, bool const viscoelastic)
{
  Medium const& ref =
    viscoelastic ? this->medium_viscoelastic_ref : this->medium_elastic_ref;

  // omega = 0 is degenerate in the constant-Q model (the velocities
  // vanish); the reference values are kept there.
  if (not this->reference_frequency_given || not(omega > 0.0)) {
    return ref;
  }

  // Kjartansson (1979): both the phase velocity and the complex velocity
  // scale as (omega / omega_0)^gamma (eq. 24), independently for the P wave
  // (with Qp) and for the S wave (with Qs).
  double const ratio = omega / (2.0 * std::numbers::pi * this->f_ref);

  size_t const n_layers = this->thickness.size();
  Medium ans = ref;
  std::vector<std::complex<double>> vs_elastic(n_layers);
  for (size_t i = 0; i < n_layers; ++i) {
    double const scale_p = std::pow(ratio, this->gamma_p[i]);
    double const scale_s = std::pow(ratio, this->gamma_s[i]);
    ans.vp[i] = ref.vp[i] * scale_p;
    ans.vs[i] = ref.vs[i] * scale_s;
    ans.mu[i] = ans.vs[i] * ans.vs[i] * this->rho[i];
    vs_elastic[i] = this->medium_elastic_ref.vs[i] * scale_s;
  }
  // The S-wave velocities of the layers change with the frequency
  // differently from each other, so do the LVLs.
  ans.LVL_indices = this->detect_LVLs(vs_elastic);

  return ans;
}

Medium
SurfaceWaveDispersion::interpolate_media(Medium const& elastic,
                                         Medium const& viscoelastic,
                                         double const fraction)
{
  Medium ans = elastic;
  for (size_t i = 0; i < ans.vs.size(); ++i) {
    ans.vp[i] =
      elastic.vp[i] * std::pow(viscoelastic.vp[i] / elastic.vp[i], fraction);
    ans.vs[i] =
      elastic.vs[i] * std::pow(viscoelastic.vs[i] / elastic.vs[i], fraction);
    ans.mu[i] = ans.vs[i] * ans.vs[i] * this->rho[i];
  }
  return ans;
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
                                    Medium const& medium)

{
  auto const target_function = [this, &medium](double phv, double omega) {
    return this->rayleigh_secular_function(phv, omega, medium);
  };

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
SurfaceWaveDispersion::nu(std::complex<double> const k,
                          double const omega,
                          int const ind,
                          Medium const& medium)
{
  std::complex<double> const k_s = omega / medium.vs[ind];
  std::complex<double> ans = std::sqrt(k * k - k_s * k_s);
  if (ans.real() < 0) {
    ans *= -1.0;
  }
  return ans;
}

std::complex<double>
SurfaceWaveDispersion::gamma(std::complex<double> const k,
                             double const omega,
                             int const ind,
                             Medium const& medium)
{
  std::complex<double> const k_p = omega / medium.vp[ind];
  std::complex<double> ans = std::sqrt(k * k - k_p * k_p);
  if (ans.real() < 0) {
    ans *= -1.0;
  }
  return ans;
}

std::complex<double>
SurfaceWaveDispersion::chi1(std::complex<double> const k,
                            double const omega,
                            int const ind,
                            Medium const& medium)
{
  std::complex<double> const k_s = omega / medium.vs[ind];
  std::complex<double> ans = 2.0 * k - k_s * k_s / k;
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
SurfaceWaveDispersion::is_LVL(int const ind, Medium const& medium)
{
  return std::find(medium.LVL_indices.begin(), medium.LVL_indices.end(), ind) !=
         medium.LVL_indices.end();
}

void
SurfaceWaveDispersion::rayleigh_downward_recursion(
  std::complex<double> const k,
  double const omega,
  Medium const& medium,
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
    std::complex<double> const mu_i = medium.mu[i];
    std::complex<double> const nu_i = this->nu(k, omega, i, medium);
    std::complex<double> const gamma_i = this->gamma(k, omega, i, medium);
    std::complex<double> const chi1_i = this->chi1(k, omega, i, medium);

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
                                                 double const omega,
                                                 Medium const& medium)
{
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>> vec_E_i;
  std::vector<Eigen::Matrix<std::complex<double>, 4, 4>> vec_E_inv_i;
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>> vec_Lambda_kernel_i;
  std::vector<Eigen::Matrix<std::complex<double>, 2, 2>> vec_R_du_i;

  this->rayleigh_downward_recursion(
    k, omega, medium, vec_E_i, vec_E_inv_i, vec_Lambda_kernel_i, vec_R_du_i);

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
                                                double const radius,
                                                Medium const& medium)
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
    residue +=
      dk * this->rayleigh_surface_response(k_n + dk, omega, medium)(1, 1);
  }
  residue /= (double)this->contour_npts;

  return k_n * residue;
}

std::vector<double>
SurfaceWaveDispersion::rayleigh_secular_function(double const phv,
                                                 double const omega,
                                                 Medium const& medium)
{
  std::complex<double> const k{ omega / phv, 0.0 };

  std::vector<double> ans;
  for (std::complex<double> const& v :
       this->rayleigh_secular_functions(k, omega, medium)) {
    ans.push_back(v.real());
    ans.push_back(v.imag());
  }
  return ans;
}

std::vector<std::complex<double>>
SurfaceWaveDispersion::rayleigh_secular_functions(std::complex<double> const k,
                                                  double const omega,
                                                  Medium const& medium)
{
  Eigen::Matrix<std::complex<double>, 2, 2> const Imat_2x2 =
    Eigen::Matrix<std::complex<double>, 2, 2>::Identity();

  std::vector<std::complex<double>> ans;

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
    k, omega, medium, vec_E_i, vec_E_inv_i, vec_Lambda_kernel_i, vec_R_du_i);

  if (n_layers >= 2) {
    Eigen::Matrix<std::complex<double>, 4, 4> E_i = vec_E_i[0];
    Eigen::Matrix<std::complex<double>, 2, 2> E21 = E_i.block<2, 2>(2, 0);
    Eigen::Matrix<std::complex<double>, 2, 2> E22 = E_i.block<2, 2>(2, 2);

    Eigen::Matrix<std::complex<double>, 2, 2> Lambda_u_1_z0 =
      vec_Lambda_kernel_i[0];

    std::complex<double> r_scf0 =
      (E21 + E22 * Lambda_u_1_z0 * vec_R_du_i[0]).determinant();
    ans.push_back(r_scf0);

    vec_R_ud_i[0] = -E21.inverse() * E22 * Lambda_u_1_z0;
  }

  if (medium.LVL_indices.empty()) {
    return ans;
  }

  // Propagate matrices (upgoing)
  // for ind = 1 to j_{max}
  // (ind = i) -> i = 1 to j_{max}

  int const max_LVL_ind = medium.LVL_indices.back();
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

    if (this->is_LVL(i, medium)) {
      std::complex<double> r_scfj =
        (Imat_2x2 - vec_R_ud_i[i] * vec_R_du_i[i]).determinant();
      ans.push_back(r_scfj);
    }
  }
  return ans;
}

void
SurfaceWaveDispersion::calculate_rayleigh_power_fraction()
{
  // Absolute value of the medium response A (Harkrider, 1964) of each mode
  // for a vertical point force and the vertical displacement at the free
  // surface, and the power of each mode relative to that of the 0th mode,
  // following the Thomson-Haskell implementation: (A_j^2 c_j) / (A_0^2 c_0).
  // Here A is evaluated as the residue of the surface response; see
  // rayleigh_medium_response().
  this->medium_responses_rayleigh.clear();
  this->psds_rayleigh.clear();

  size_t n_unresolved = 0;
  size_t n_invalid_0th = 0;

  for (size_t i = 0; i < this->phvs_rayleigh.size(); ++i) {
    std::vector<double> const& phvs_at_freq = this->phvs_rayleigh[i];
    std::vector<double> psds_at_freq(phvs_at_freq.size(), std::nan(""));
    std::vector<double> medium_responses_at_freq(phvs_at_freq.size(),
                                                 std::nan(""));

    // f = 0: the medium response is not defined.
    if (this->freq[i] == 0.0) {
      this->medium_responses_rayleigh.push_back(medium_responses_at_freq);
      this->psds_rayleigh.push_back(psds_at_freq);
      continue;
    }

    double const omega = this->freq[i] * 2.0 * std::numbers::pi;
    Medium const medium = this->medium_at(omega, false);

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
        this->rayleigh_medium_response(phv, omega, radius, medium);
      double const power = std::norm(medium_response) * phv;
      medium_responses_at_freq[j] = std::abs(medium_response);

      if (j == 0) {
        power_0th = power;
      }
      psds_at_freq[j] = power / power_0th;
    }

    if (not phvs_at_freq.empty() &&
        not(std::isfinite(power_0th) && power_0th > 0.0)) {
      ++n_invalid_0th;
    }
    this->medium_responses_rayleigh.push_back(medium_responses_at_freq);
    this->psds_rayleigh.push_back(psds_at_freq);
  }

  if (n_unresolved > 0) {
    std::cerr << "[Rayl medium response] " << n_unresolved
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
SurfaceWaveDispersion::predicted_num_of_modes(double phv,
                                              double omega,
                                              Medium const& medium)
{
  double ans = 0.0;
  for (size_t i = 0; i < this->thickness.size() - 2; ++i) {
    std::complex<double> gam_p = std::sqrt(
      std::complex<double>{ std::pow(phv / medium.vp[i].real(), 2) - 1.0, 0 });

    std::complex<double> gam_s = std::sqrt(
      std::complex<double>{ std::pow(phv / medium.vs[i].real(), 2) - 1.0, 0 });

    ans += ((gam_s + gam_p) * this->thickness[i]).real();
  }
  // ans = omega / (2 * pi) * 2 * (sum (gamma + nu)*h) / phv
  ans *= omega / (std::numbers::pi * phv);

  return ans;
}

bool
SurfaceWaveDispersion::equals_vs_of_a_layer(double phv, Medium const& medium)
{
  bool ans = false;

  for (size_t i = 0; i < medium.vs.size(); ++i) {
    if (std::abs(phv - medium.vs[i].real()) < 1e-3) {
      ans = true;
    }
  }

  return ans;
}

bool
SurfaceWaveDispersion::equals_vp_of_a_layer(double phv, Medium const& medium)
{
  bool ans = false;

  for (size_t i = 0; i < medium.vp.size(); ++i) {
    if (std::abs(phv - medium.vp[i].real()) < 1e-3) {
      ans = true;
    }
  }

  return ans;
}

void
SurfaceWaveDispersion::calculate_rayleigh_phase_velocities()
{
  this->phvs_rayleigh.clear();
  this->keys_rayleigh.clear();

  for (size_t i = 0; i < this->freq.size(); ++i) {
    std::vector<double> phvs_rayleigh_at_freq;
    std::vector<int> keys_rayleigh_at_freq;

    double omega = freq[i] * 2.0 * std::numbers::pi;

    // The root search uses the phase velocities of the elastic medium at
    // this frequency.
    Medium const medium = this->medium_at(omega, false);
    double const vs_half_space = medium.vs.back().real();

    if (this->freq[i] == 0.0) {
      // f = 0 is degenerate; this row is extrapolated afterwards.
      this->phvs_rayleigh.push_back(std::vector<double>{ vs_half_space });
      this->keys_rayleigh.push_back(std::vector<int>{ 0 });
      continue;
    }

    // Predicte number of modes using Fan et al. (2007), published by Chinese
    // Journal of Geophysics. Note that it depends on frequency.
    int const num_of_grid_pts =
      (this->predicted_num_of_modes(vs_half_space, omega, medium) + 1) *
      this->ngrid;

    double min_vs = medium.vs[0].real();
    for (size_t j = 1; j < medium.vs.size(); ++j) {
      min_vs = std::min(min_vs, medium.vs[j].real());
    }
    double const dc = (vs_half_space - 0.8 * min_vs) / (double)num_of_grid_pts;

    // Prepare grid points for grid search
    std::vector<double> grid_pts;
    grid_pts.resize(num_of_grid_pts);
    for (int j = 0; j < num_of_grid_pts; ++j) {
      grid_pts[j] =
        0.8 * min_vs + (double)(j + 1) * dc + this->fluctuation_for_root_search;
    }
    if (i > 0 && this->freq[i - 1] > 0.0) {
      for (size_t j = 0; j < phvs_rayleigh[i - 1].size(); ++j) {
        double center = phvs_rayleigh[i - 1][j];

        //  Following an idea presented by Fan et al. (2007).
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
    std::vector<double> bf = rayleigh_secular_function(phv_bf, omega, medium);
    for (size_t j = 1; j < grid_pts.size(); ++j) {
      double phv_af = grid_pts[j];
      std::vector<double> af = rayleigh_secular_function(phv_af, omega, medium);

      // Check kth secular function
      for (size_t k = 0; k < af.size(); ++k) {
        if (bf[k] * af[k] < 0 && std::abs(bf[k] * af[k]) > 1e-8) {
          double candidate = bisection_sf(phv_bf, phv_af, omega, k, medium);
          if (not std::isnan(candidate) && not std::isinf(candidate) &&
              not this->equals_vs_of_a_layer(candidate, medium) &&
              not this->equals_vp_of_a_layer(candidate, medium)) {
            phvs_rayleigh_at_freq.push_back(candidate);
            keys_rayleigh_at_freq.push_back(k);
            break;
          }
        }
      }
      phv_bf = phv_af;
      bf = std::move(af);
    }
    this->phvs_rayleigh.push_back(phvs_rayleigh_at_freq);
    this->keys_rayleigh.push_back(keys_rayleigh_at_freq);
  }

  if (this->freq[0] == 0.0) {
    if (phvs_rayleigh[1].size() > 0 and phvs_rayleigh[2].size() > 0) {
      this->phvs_rayleigh[0][0] =
        this->phvs_rayleigh[1][0] * 2 - this->phvs_rayleigh[2][0];
    } else {
      std::cerr << "Something went wrong in low frequency region." << std::endl;
      std::cerr << "Increasing ngrid might resolve the problem." << std::endl;
      this->phvs_rayleigh[0][0] = this->phvs_rayleigh[1][0];
    }
  }
  inspect_nmodes_rayleigh();
}

bool
SurfaceWaveDispersion::newton_root(std::complex<double>& k,
                                   double const omega,
                                   int const sf_index,
                                   Medium const& medium)
{
  // The secular functions are analytic in the complex wavenumber (away from
  // the branch cuts of nu and gamma), so that the derivative is taken along
  // the real axis.
  auto const secular_function = [&](std::complex<double> const k_val) {
    return this->rayleigh_secular_functions(k_val, omega, medium)[sf_index];
  };

  for (int itr = 0; itr < this->newton_max_itr; ++itr) {
    double const h = this->newton_rel_diff_step * std::abs(k);
    std::complex<double> const f = secular_function(k);
    std::complex<double> const derivative = (secular_function(k + h) - f) / h;

    std::complex<double> dk = -f / derivative;
    if (not(std::isfinite(dk.real()) && std::isfinite(dk.imag()))) {
      return false;
    }

    double const limit = this->newton_rel_step_limit * std::abs(k);
    if (std::abs(dk) > limit) {
      dk *= limit / std::abs(dk);
    }
    k += dk;

    if (std::abs(dk) < this->newton_rel_tol * std::abs(k)) {
      return true;
    }
  }
  return false;
}

bool
SurfaceWaveDispersion::is_acceptable_root(std::complex<double> const k,
                                          double const k_elastic)
{
  double const alpha = -k.imag();
  // A tiny negative alpha is a rounding error of alpha = 0 (elastic medium).
  return k.real() > (1.0 - this->wavenumber_rel_range) * k_elastic &&
         k.real() < (1.0 + this->wavenumber_rel_range) * k_elastic &&
         alpha > -1e-9 * k_elastic && alpha < this->alpha_max;
}

std::complex<double>
SurfaceWaveDispersion::track_root_from_elastic(std::complex<double> k,
                                               double const omega,
                                               int const sf_index,
                                               Medium const& elastic,
                                               Medium const& viscoelastic)
{
  // Starting directly from alpha = 0 with the full attenuation can end up at
  // another root or at a point that is not a root. The attenuation is
  // therefore switched on gradually and the root is followed. The solution
  // before the last step (full attenuation) is returned.
  for (int step = 1; step < this->tracking_num_steps; ++step) {
    double const fraction = (double)step / (double)this->tracking_num_steps;
    Medium const medium =
      this->interpolate_media(elastic, viscoelastic, fraction);

    std::complex<double> k_next = k;
    if (this->newton_root(k_next, omega, sf_index, medium)) {
      k = k_next;
    }
  }
  return k;
}

void
SurfaceWaveDispersion::interpolate_missing_modes(
  std::vector<std::vector<double>>& table)
{
  size_t const n_freq = table.size();
  size_t const n_modes = this->nmodes_rayleigh;

  // Index of the frequency at which each mode cuts in.
  std::vector<size_t> cut_in(n_modes, n_freq);
  for (size_t j = 0; j < n_modes; ++j) {
    for (size_t i = 0; i < n_freq; ++i) {
      if (not std::isnan(table[i][j])) {
        cut_in[j] = i;
        break;
      }
    }
  }

  // A frequency at which a mode is missing after its cut-in is invalid for
  // all the modes, since a missing root shifts the higher modes by one
  // column.
  std::vector<bool> valid(n_freq, true);
  for (size_t j = 0; j < n_modes; ++j) {
    for (size_t i = cut_in[j]; i < n_freq; ++i) {
      if (std::isnan(table[i][j])) {
        valid[i] = false;
      }
    }
  }

  for (size_t j = 0; j < n_modes; ++j) {
    std::vector<double> xs, ys;
    for (size_t i = cut_in[j]; i < n_freq; ++i) {
      if (valid[i]) {
        xs.push_back(this->freq[i]);
        ys.push_back(table[i][j]);
      }
    }
    if (xs.size() < 2) {
      continue;
    }

    for (size_t i = cut_in[j]; i < n_freq; ++i) {
      if (valid[i]) {
        continue;
      }
      // Linear interpolation (extrapolation outside the valid range).
      size_t seg =
        std::upper_bound(xs.begin(), xs.end(), this->freq[i]) - xs.begin();
      seg = std::clamp(seg, (size_t)1, xs.size() - 1) - 1;
      double const slope = (ys[seg + 1] - ys[seg]) / (xs[seg + 1] - xs[seg]);
      table[i][j] = slope * (this->freq[i] - xs[seg]) + ys[seg];
    }
  }
}

void
SurfaceWaveDispersion::calculate_rayleigh_attenuation()
{
  // The complex roots k = omega / phv - i * alpha of the secular functions
  // of the viscoelastic medium are searched for by the Newton method. The
  // initial guess of a mode is the root of the elastic medium and the
  // attenuation coefficient of the same mode at the previous frequency.
  double const nan = std::nan("");
  size_t const n_freq = this->freq.size();
  size_t const n_modes = this->nmodes_rayleigh;

  std::vector<std::vector<double>> phvs(n_freq,
                                        std::vector<double>(n_modes, nan));
  std::vector<std::vector<int>> keys(n_freq, std::vector<int>(n_modes, 0));
  for (size_t i = 0; i < n_freq; ++i) {
    for (size_t j = 0; j < this->phvs_rayleigh[i].size(); ++j) {
      phvs[i][j] = this->phvs_rayleigh[i][j];
      if (j < this->keys_rayleigh[i].size()) {
        keys[i][j] = this->keys_rayleigh[i][j];
      }
    }
  }
  this->interpolate_missing_modes(phvs);

  std::vector<std::vector<double>> attens(n_freq,
                                          std::vector<double>(n_modes, nan));

  size_t n_unresolved = 0;

  for (size_t i = 0; i < n_freq; ++i) {
    double const omega = this->freq[i] * 2.0 * std::numbers::pi;
    if (not(omega > 0.0)) {
      // f = 0 is degenerate; this row is extrapolated afterwards.
      continue;
    }

    Medium const elastic = this->medium_at(omega, false);
    Medium const viscoelastic = this->medium_at(omega, true);
    int const n_sf = 1 + (int)viscoelastic.LVL_indices.size();

    std::vector<bool> exists(n_modes, false);
    std::vector<bool> solved(n_modes, false);
    std::vector<double> k_elastic(n_modes, nan);
    std::vector<std::complex<double>> k_start(n_modes);
    std::vector<std::complex<double>> k_root(n_modes);

    for (size_t j = 0; j < n_modes; ++j) {
      if (std::isnan(phvs[i][j])) {
        continue;
      }
      exists[j] = true;
      k_elastic[j] = omega / phvs[i][j];

      double alpha_init = 0.0;
      if (i > 0 && not std::isnan(attens[i - 1][j])) {
        alpha_init = attens[i - 1][j];
      }

      // The secular function with which the root was found.
      int sf_index = keys[i][j] / 2;
      if (sf_index >= n_sf) {
        sf_index = 0;
      }

      std::complex<double> k{ k_elastic[j], -alpha_init };
      if (i == 0) {
        // First row of a frequency range that does not start from 0 Hz:
        // there is no previous frequency to take the attenuation from.
        k = this->track_root_from_elastic(
          k, omega, sf_index, elastic, viscoelastic);
      }
      k_start[j] = k;

      if (this->newton_root(k, omega, sf_index, viscoelastic) &&
          this->is_acceptable_root(k, k_elastic[j])) {
        solved[j] = true;
        k_root[j] = k;
      }
    }

    // Where two modes cross (which they can with attenuation), one of them
    // is confined around an LVL and the other one is not. The root of the
    // former then almost cancels with a pole in the secular function of the
    // free surface, and that of the latter in the secular function of the
    // LVL, so that each secular function leads to only one of the two roots
    // from both initial guesses. A mode that converged to the same root as
    // another one, or that was not solved, is therefore solved again with
    // every secular function, and a root that differs from those of all the
    // other modes is taken.
    for (size_t j = 0; j < n_modes; ++j) {
      if (not exists[j]) {
        continue;
      }

      std::vector<size_t> start_modes{ j };
      if (solved[j]) {
        for (size_t jp = 0; jp < j; ++jp) {
          if (solved[jp] && std::abs(k_root[j] - k_root[jp]) <
                              this->coincidence_rel_tol * std::abs(k_root[j])) {
            start_modes.push_back(jp);
          }
        }
        if (start_modes.size() == 1) {
          continue;
        }
      }

      bool found = false;
      double best_distance = std::numeric_limits<double>::infinity();
      std::complex<double> best_k;
      int best_sf_index = 0;

      for (int sf_index = 0; sf_index < n_sf; ++sf_index) {
        for (size_t const jp : start_modes) {
          std::complex<double> k = k_start[jp];
          if (not this->newton_root(k, omega, sf_index, viscoelastic) ||
              not this->is_acceptable_root(k, k_elastic[jp])) {
            continue;
          }

          bool coincides = false;
          for (size_t jq = 0; jq < n_modes; ++jq) {
            if (jq != j && solved[jq] &&
                std::abs(k - k_root[jq]) <
                  this->coincidence_rel_tol * std::abs(k)) {
              coincides = true;
            }
          }
          if (coincides) {
            continue;
          }

          // Distance to the root estimated by a Newton step, relative to
          // the wavenumber.
          double const h = this->newton_rel_diff_step * std::abs(k);
          std::complex<double> const derivative =
            (this->rayleigh_secular_functions(
               k + h, omega, viscoelastic)[sf_index] -
             this->rayleigh_secular_functions(
               k - h, omega, viscoelastic)[sf_index]) /
            (2.0 * h);
          double const distance = std::abs(this->rayleigh_secular_functions(
                                             k, omega, viscoelastic)[sf_index] /
                                           derivative) /
                                  std::abs(k);
          if (distance < this->coincidence_rel_tol &&
              distance < best_distance) {
            found = true;
            best_distance = distance;
            best_k = k;
            best_sf_index = sf_index;
          }
        }
      }

      if (found) {
        solved[j] = true;
        k_root[j] = best_k;
        keys[i][j] = 2 * best_sf_index;
      } else {
        solved[j] = false;
        ++n_unresolved;
      }
    }

    for (size_t j = 0; j < n_modes; ++j) {
      if (not exists[j]) {
        continue;
      }
      if (solved[j]) {
        phvs[i][j] = omega / k_root[j].real();
        attens[i][j] = std::max(-k_root[j].imag(), 0.0);
      } else {
        phvs[i][j] = nan;
        attens[i][j] = nan;
      }
    }

    // Modes are numbered in ascending order of the phase velocity at each
    // frequency, as in the elastic case. With attenuation two modes can
    // cross, which swaps them here.
    std::vector<size_t> order(n_modes);
    for (size_t j = 0; j < n_modes; ++j) {
      order[j] = j;
    }
    std::stable_sort(
      order.begin(), order.end(), [&](size_t const a, size_t const b) {
        return not std::isnan(phvs[i][a]) &&
               (std::isnan(phvs[i][b]) || phvs[i][a] < phvs[i][b]);
      });
    std::vector<double> const phvs_at_freq = phvs[i];
    std::vector<double> const attens_at_freq = attens[i];
    std::vector<int> const keys_at_freq = keys[i];
    for (size_t j = 0; j < n_modes; ++j) {
      phvs[i][j] = phvs_at_freq[order[j]];
      attens[i][j] = attens_at_freq[order[j]];
      keys[i][j] = keys_at_freq[order[j]];
    }
  }

  if (this->freq[0] == 0.0 && n_freq > 2) {
    // f = 0 is degenerate; extrapolate from the next two points.
    attens[0][0] = std::max(attens[1][0] * 2 - attens[2][0], 0.0);
  }

  if (n_unresolved > 0) {
    std::cerr << "[Rayl attenuation] " << n_unresolved
              << " mode(s) did not converge to a root in the accepted range "
                 "or could not be separated from another mode (written as "
                 "nan)."
              << std::endl;
  }

  this->phvs_rayleigh = phvs;
  this->keys_rayleigh = keys;
  this->attens_rayleigh = attens;
}
