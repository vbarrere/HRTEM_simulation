# HRTEM_simulation

Multislice HRTEM image simulation of MD nanoparticles (AgCo) on amorphous-carbon substrates, parallelised with MPI.
Settings and paths are in `run.sh`; `make` builds `main`.

## Requirements

- An MPI Fortran compiler (`mpifort`), FFTW3 (double precision) and OpenBLAS. The library paths are set in `LDLIBS`
  in the `Makefile` (`/lib/x86_64-linux-gnu/libfftw3.so.3` and `/lib/x86_64-linux-gnu/libopenblas.so.0`, from the
  Ubuntu packages `libfftw3-double3` and `libopenblas0`).
- `run.sh` sets `OPENBLAS_NUM_THREADS=1`. Each MPI rank must use a single BLAS thread; otherwise every rank starts one
  thread per core and the machine is oversubscribed.

## Output

A run writes a single file, `output_data` in `run.sh` (`dataset.npy`). It is a NumPy file holding one record per
image, and each record contains the descriptors and the image of that simulation. Both are written by one binary
`write` statement (`save_data_row` in `utils_io.f90`), so a descriptor row cannot be paired with the wrong image.

Fields of a record, in file order:

| Fields | Type |
|---|---|
| `id_sim n_atoms n_steps nat1 nat2 nat1_out nat2_out nat1_in nat2_in` | int32 |
| `initial_temperature epot_total composition gyration_radius d_com coreshell_index` | float64 |
| `phi theta psi position_scale` (rotation angles in rad, scale factor) | float64 |
| `defocus A1_re A1_im B2_re B2_im A2_re A2_im C3 S3_re S3_im A3_re A3_im B4_re B4_im D4_re D4_im A4_re A4_im` (nm) | float64 |
| `image` (ny × nx, first index = y, min-max scaled to 0-255) | uint8 |

`defocus` is C1. `d_com` is NaN when the cluster contains a single species. A record takes 260 + n² bytes: 64.3 KB
at 256 px. The field names and their order are set in `int_fields` and `real_fields` (`utils_io.f90`).

```python
import numpy as np
data = np.load("dataset.npy", mmap_mode="r")      # memory-mapped, nothing is read until used
data["image"][0]                                  # first image, shape (256, 256), dtype uint8
data[0]["phi"], data[0]["defocus"]                # descriptors of the same record
data.dtype.names                                  # all field names

import pandas as pd                               # descriptors as a table
descriptors = pd.DataFrame({name: data[name] for name in data.dtype.names if name != "image"})
```

Reading one descriptor for all images touches every record of the file, which is slow for a large dataset. If the
descriptors are needed often, extract them once (as above) and save that table; `id_sim` links it back to the
records.

`images_to_png.py` writes one `{id_sim}.png` per image if PNG files are needed. For training, reading the
memory-mapped array directly avoids creating one file per image.

## Debye-Waller factors

The isotropic Debye-Waller factors are set in `constants.f90`, in nm² (1 nm² = 100 Å²):

```fortran
double precision, parameter :: dbf_ag = 0.0d0, dbf_co = 0.0d0, dbf_c = 0.001442d0 ! nm²
```

Convention: B = 8π² ⟨u_x²⟩, where ⟨u_x²⟩ is the mean-square displacement along one axis. The elastic form factor of each
atom is damped by exp(-B g²/4), with g in nm⁻¹ (`scattering_factor` in `slc.f90`). B also sets the absorptive part of the
form factor (thermal diffuse scattering); with B = 0 there is neither damping nor absorption.

The thermal motion of each atom must be counted once, either in its coordinates or in B. The particle and the substrate
come from different kinds of MD output, so they are treated differently:

| Atoms | Coordinates | B |
|---|---|---|
| Ag, Co (particle) | one instantaneous MD snapshot | 0 (the displacements are in the coordinates) |
| C (a-C substrate) | averaged over 500 fs at 300 K | 0.1442 Å² = 0.001442 nm² (diamond at 300 K) |

Reference for the tabulated values below: L.-M. Peng, G. Ren, S. L. Dudarev and M. J. Whelan, "Debye-Waller factors and
absorptive scattering factors of elemental crystals", Acta Cryst. A52, 456-470 (1996), doi:10.1107/S010876739600089X,
supplementary Table 1 (B of 44 elemental crystals from 1 K upwards, computed from measured phonon densities of states).

### Particle (Ag, Co): B = 0

Each particle is a single instantaneous MD snapshot, so its coordinates contain the thermal displacements at the MD
temperature: each image is one frozen-phonon configuration, and a B on top would count the thermal motion twice.

- This holds if the MD ends at the temperature the images represent (300 K). Classically the mean-square displacement is
  proportional to T: a snapshot taken at 600 K carries twice the displacement variance of one at 300 K.
- Classical MD has no zero-point motion. For Ag and Co at 300 K this changes B by about +2 % and +6 % (Debye model).
- An experimental image averages over the configurations visited during the exposure. In the linear (weak-phase)
  approximation this damps each spatial frequency of the averaged structure by exp(-B g²/4), which is the same as a
  Gaussian image blur of rms √(B/8π²) per axis: about 0.010 nm for Ag (B = 0.75 Å² at 300 K, Peng et al.) and 0.006 nm
  for Co (B ≈ 0.3 Å², Debye model). One snapshot is not averaged. In a crystal, the sum over the many atoms of a lattice
  plane averages the displacements, so a snapshot already shows damped fringes plus a weak diffuse background. In a
  disordered particle there is no such averaging: the high-frequency contrast keeps its full strength, and a fraction
  1 - exp(-B g²/2) of it is set by the displacements instead of the structure (Ag: 6.5 % at 0.236 nm, 24 % at the
  0.117 nm band limit). The snapshot then has up to exp(B g²/4) more contrast amplitude than the averaged image: +3 % at
  0.236 nm and +15 % at 0.117 nm for Ag, +2 % and +6 % for Co.
- This difference lies within the image blur sampled in `read_input` (`vib1`, `vib2`: Gaussian rms 0.005-0.04 nm per
  direction, `apply_vibration` in `wavimg.f90`), which damps the band limit by factors between 0.97 and 0.10. Adding the
  thermal blur in quadrature only moves the lower end of that range, from 0.005 nm to 0.011 nm for Ag.
- The snapshot also carries the larger vibrations of surface atoms, which a single B per element cannot represent.

### Substrate (C): B = 0.144 Å²

The value is B of diamond at 300 K in Peng et al. (0.1442 Å²; 0.1435 Å² at 293 K).

- For carbon, B at room temperature is mostly zero-point motion: diamond has B = 0.1286 Å² at 1 K. Classical MD omits
  it, so the per-atom `dwf` column of the slab files (classical, slab means 0.14-0.22 Å²) is not used; the substrate
  reader ignores it.
- Diamond is a lower bound for a-C, which is mostly sp²: the classical slab values are 2-4 times the classical value of
  diamond (about 0.06 Å² at 300 K, Debye model), and a Debye-model zero-point correction of the slab values gives
  0.21-0.30 Å².
- The image is not sensitive to this choice: at the 256 px band limit (8.53 nm⁻¹, 0.117 nm), exp(-B g²/4) is 0.974 for
  0.144 Å² and 0.947 for 0.3 Å².

#### Substrate positions are time-averaged (verified)

The slabs are made with LAMMPS (Tersoff carbon, melt-quench to 300 K; `in.substrate` / `in.substrate.old-lammps` and
`run.sh` of the `carbon_substrate` generator). Their positions are time averages, so applying B to them does not count
the thermal motion twice:

- During the last 500 fs at 300 K, `fix ave/atom 1 500 500` averages the unwrapped coordinates and their squares at
  every step. The file holds the averaged positions ⟨r⟩ and, per atom, B = 8π²/3 (⟨r²⟩ - ⟨r⟩²), i.e. with the factor
  1/3 that turns the 3D mean-square displacement into ⟨u_x²⟩.
- Checked on a 5 x 5 x 2 nm test run that saved a frame every 10 fs through the averaging window: the positions in
  the output agree with an independent average of the raw frames to 0.001 Å (rms) and differ from the instantaneous
  last frame by 0.092 Å (rms). Mean B: 0.1986 Å² in the file, 0.1987 Å² recomputed, per-atom correlation 1.000.
- Checked on 21 production slabs: each file holds the averaged frame (not the empty step-0 frame LAMMPS also writes),
  the `dwf_mean` in the header matches the column, z lies in [0, thickness], and x, y go up to 0.18 Å outside
  [0, 100) (unwrapped averages; the HRTEM code wraps them when it applies the random lateral shift).
- So the slab means of 0.15-0.20 Å² are correct classical values, without zero-point motion.
- Minor: 0.02-0.35 % of the atoms hop between two sites during the window (B > 1 Å²) and their averaged position lies
  between the two; the 500 fs window does not fully average the slowest in-plane vibrations (period about 1 ps at
  10 nm wavelength). Both are negligible at 0.12 nm resolution.

### Former values for Ag and Co (debye.py)

The code used dbf_ag = 0.019 nm² and dbf_co = 0.015 nm² until October 2026. They come from `debye.py` (repository
`vbarrere/Debye-Waller_factor`), which computes 8π² ⟨|r - ⟨r⟩|²⟩ from a series of MD frames. ⟨|r - ⟨r⟩|²⟩ sums the
three axes, so this is 3 B: B = 8π² ⟨|r - ⟨r⟩|²⟩ / 3. Divided by 3, Ag gives 0.63 Å² and Co 0.5 Å² (bulk Ag at 300 K:
0.75 Å² in Peng et al.; Co is not tabulated there). Such values apply only to time-averaged particle coordinates. (Minor:
the loop over frames starts at the second frame but divides by the total number of frames.)

### Form factor cache: one B per element

The form factors (Weickenmeier-Kohl, with absorptive part) are computed once per run for each distinct (Z, B) pair and
reused for every image (`unique_scattering_factors` and `compute_form_factor` in `slc.f90`). This is valid because the
grid, the voltage and B are fixed for the whole run. At 256 px the computation takes 1-2 s per pair and per rank,
once per run instead of once per image.

Consequences:

- **B must be a single value per element.** The cache has `n_types_max` = 5 slots. A per-slab or per-atom B would fill
  it and stop the run with `more than n_types_max distinct (Z, B) pairs`. A per-atom B would need a different scheme
  (for example a few B bins per element).
- To change a value, edit `constants.f90` and rebuild. Nothing else depends on it.
