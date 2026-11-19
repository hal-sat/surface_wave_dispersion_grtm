# surface_wave_dispersion_grtm

Phase velocities and relative modal powers of Rayleigh waves in a layered
elastic half-space, computed by the generalized reflection/transmission (R/T)
coefficient method. The theory is summarized in `note/main.pdf`.

## Requirements

- g++ (C++20) and Eigen 3 (expected in `/usr/include/eigen3`)
- Python 3 with NumPy

## Build

```bash
make
```

## Usage

```bash
python3 run.py -if demo/OrgM1.profile
```

The results are written to `<name>_results_grtm/` next to the input file
(`-od` changes the directory). Other options: `-nfp` (number of frequency
points), `-ubf` (upper bound of frequency in Hz) and `-ngp` (number of grid
points for the root search).

## Input

One line per layer, from the top; the last line is the half-space.

```
# Thickness (m), Vp (m/s), Vs (m/s), rho (t/m^3)
40, 800, 300, 1.5
infty, 1500, 700, 1.8
```

## Output

- `Rayl_phvs.csv`: frequency (Hz) followed by the phase velocity of each mode
- `Rayl_psds.csv`: frequency (Hz) followed by the power of each mode relative
  to that of the mode with the lowest phase velocity
- `employed_params.txt`: parameters used for the computation

Modes that do not exist at a frequency are written as `nan`.

## License

MIT License. See `LICENSE`.
