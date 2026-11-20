// SPDX-FileCopyrightText: Copyright 2024, 2026 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "surface_wave_dispersion.hpp"

// Usage: main <profile> <result directory> <nf> <fmax> <ngrid>
//             [<fmin> <f0> [<medium>]]
// <f0> is the reference frequency (Hz) of the velocities in the profile or
// "none" (frequency-independent complex velocities).
// <medium> is "auto" (viscoelastic if Qp and Qs are given in the profile) or
// "elastic" (Qp and Qs are ignored).
int
main(int argc, char** argv)
{
  if (argc != 6 && argc != 8 && argc != 9) {
    std::cerr << "Error: Invalid argument vector." << std::endl;
    std::cerr << "Usage: " << argv[0]
              << " <profile> <result directory> <nf> <fmax> <ngrid> "
                 "[<fmin> <f0> [<medium>]]"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string inf = argv[1];
  std::string resdir = argv[2];
  int nf = std::stod(argv[3]);
  double fmax = std::stod(argv[4]);
  int ngrid = std::stoi(argv[5]);

  double fmin = 0.0;
  bool reference_frequency_given = false;
  double f_ref = 0.0;
  bool ignore_attenuation = false;
  if (argc >= 8) {
    fmin = std::stod(argv[6]);
    if (std::string{ argv[7] } != "none") {
      reference_frequency_given = true;
      f_ref = std::stod(argv[7]);
    }
  }
  if (argc == 9) {
    std::string const medium = argv[8];
    if (medium == "elastic") {
      ignore_attenuation = true;
    } else if (medium != "auto") {
      std::cerr << "Error: <medium> must be either \"auto\" or \"elastic\"."
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (nf < 2) {
    std::cerr << "Error: nf must be 2 or larger." << std::endl;
    return EXIT_FAILURE;
  }
  if (fmin < 0.0 || not(fmax > fmin)) {
    std::cerr << "Error: 0 <= fmin < fmax must hold." << std::endl;
    return EXIT_FAILURE;
  }
  if (ngrid < 1) {
    std::cerr << "Error: ngrid must be positive." << std::endl;
    return EXIT_FAILURE;
  }
  if (reference_frequency_given && not(f_ref > 0.0)) {
    std::cerr << "Error: The reference frequency must be positive."
              << std::endl;
    return EXIT_FAILURE;
  }

  // ---

  SurfaceWaveDispersion swd{ nf, fmin, fmax };
  swd.set_number_of_grids(ngrid);
  swd.set_ignore_attenuation(ignore_attenuation);

  swd.dump_basic_parameters();

  swd.read_profile(inf);

  if (reference_frequency_given) {
    if (not swd.has_attenuation()) {
      std::cerr << "Error: The reference frequency is meaningful only for a "
                   "viscoelastic medium (a profile with Qp and Qs that are "
                   "not ignored)."
                << std::endl;
      return EXIT_FAILURE;
    }
    swd.set_reference_frequency(f_ref);
  }

  swd.dump_ground_structure();

  swd.calculate_rayleigh_phase_velocities();
  if (swd.has_attenuation()) {
    // Viscoelastic medium: phase velocities and attenuation coefficients.
    // The relative power is not computed.
    swd.calculate_rayleigh_attenuation();
    swd.dump_number_of_modes_rayleigh();
    swd.write_rayleigh_results(resdir + "Rayl_phvs.csv");
    swd.write_rayleigh_attenuation(resdir + "Rayl_attens.csv");
  } else {
    // Elastic medium: phase velocities and relative powers.
    swd.dump_number_of_modes_rayleigh();
    swd.write_rayleigh_results(resdir + "Rayl_phvs.csv");
    swd.calculate_rayleigh_power_fraction();
    swd.write_rayleigh_power_fraction(resdir + "Rayl_psds.csv");
  }
  swd.dump_ground_structure();
  swd.dump_LVL_indices();
}
