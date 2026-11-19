// SPDX-FileCopyrightText: Copyright 2024, 2026 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "surface_wave_dispersion.hpp"

int
main(int argc, char** argv)
{
  if (argc != 6) {
    perror("Error: Invalid argument vector.");
    exit(errno);
  }

  std::string inf = argv[1];
  std::string resdir = argv[2];
  int nf = std::stod(argv[3]);
  double fmax = std::stod(argv[4]);
  int ngrid = std::stoi(argv[5]);

  double df = fmax / (double)(nf + 1);

  // ---

  SurfaceWaveDispersion swd{ nf, fmax, df };
  swd.set_number_of_grids(ngrid);

  swd.dump_basic_parameters();

  swd.read_profile(inf);

  swd.dump_ground_structure();

  swd.calculate_rayleigh_phase_velocities();
  swd.dump_number_of_modes_rayleigh();
  swd.write_rayleigh_results(resdir + "Rayl_phvs.csv");
  swd.calculate_rayleigh_power_fraction();
  swd.write_rayleigh_power_fraction(resdir + "Rayl_psds.csv");
  swd.dump_ground_structure();
  swd.dump_LVL_indices();
}
