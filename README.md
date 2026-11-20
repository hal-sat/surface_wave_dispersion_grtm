# surface_wave_dispersion_grtm

Phase velocities of Rayleigh waves in a layered half-space, computed by the
generalized reflection/transmission (R/T) coefficient method. For an elastic
medium the medium responses and the relative modal powers are also computed,
and for a viscoelastic medium (a profile with Qp and Qs) the attenuation
coefficients are computed instead. The theory is summarized in `note/main.pdf`.

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

`run.py` can be called from any directory; the input file and `-od` are
relative to the current directory.

The results are written to `<name>_results_grtm/` next to the input file
(`-od` changes the directory). Other options:

- `-nfp`: number of frequency intervals `nf`; the frequencies are
  `fmin + (fmax - fmin) / nf * i` for `i = 0` to `nf`
- `-lbf`, `-ubf`: lower and upper bounds of frequency `fmin` and `fmax` in Hz
- `-ngp`: number of grid points per predicted mode for the root search
- `-nmr`: do not compute the medium responses and the relative powers
  (elastic medium only), which saves time when only the phase velocities are
  needed
- `-f0`: reference frequency in Hz (viscoelastic medium only; see below)
- `-el`: ignore Qp and Qs of the profile and treat the medium as elastic

## Input

One line per layer, from the top; the last line is the half-space.

```
# Thickness (m), Vp (m/s), Vs (m/s), rho (t/m^3)
40, 800, 300, 1.5
infty, 1500, 700, 1.8
```

For a viscoelastic medium, add Qp and Qs to every layer (Q <= 0 means that
the layer is elastic):

```
# Thickness (m), Vp (m/s), Vs (m/s), rho (t/m^3), Qp, Qs
40, 800, 300, 1.5, 10, 10
infty, 1500, 700, 1.8, 30, 30
```

```bash
python3 run.py -if demo/OrgM1_viscoelastic.profile
```

With `-el`, Qp and Qs are ignored and the medium is treated as elastic, i.e.,
the result is that of the profile without Qp and Qs.

By default the complex velocities are frequency independent. With `-f0`, Vp
and Vs of the profile are regarded as the phase velocities at that frequency
and the causal constant-Q model of Kjartansson (1979) is used, so that the
velocities depend on the frequency.

## Output

- `Rayl_phvs.csv`: frequency (Hz) followed by the phase velocity of each mode
- `Rayl_medium_responses.csv` (elastic medium only): frequency (Hz) followed
  by the absolute value of the medium response `A` of each mode for a vertical
  point force and the vertical displacement at the free surface. Its unit
  follows from those of the input: 1/(kPa m) for m, m/s and t/m^3. The power
  of a mode is proportional to `c A^2`, where `c` is its phase velocity
- `Rayl_psds.csv` (elastic medium only): frequency (Hz) followed by the power
  of each mode relative to that of the mode with the lowest phase velocity,
  i.e. the ratio of `c A^2`
- `Rayl_attens.csv` (viscoelastic medium only): frequency (Hz) followed by the
  attenuation coefficient (1/m) of each mode
- `employed_params.txt`: parameters used for the computation

Modes are numbered in ascending order of the phase velocity at each frequency.
Modes that do not exist at a frequency, or whose root was not found, are
written as `nan`.

## License

MIT License. See `LICENSE`.
