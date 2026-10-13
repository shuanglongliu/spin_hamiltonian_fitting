# anisotropic_exchange

Extract full 3×3 exchange coupling matrices **J** from noncollinear DFT energies.

The fit is done by `main.jl` (Julia). The Python scripts are helpers for preparing its input and
post-processing its output.

---

## Spin Hamiltonian

$$H = \text{factor} \times \sum_{\langle ij \rangle} \mathbf{S}_i^{\mathsf T}\, \mathbf{J}_{ij}\, \mathbf{S}_j$$

Each **J**<sub>*ij*</sub> is expressed in a local reference frame attached to the pair (*i*, *j*).
`factor` is `−2`. Energies are read in eV and the *J* components are reported in meV.

For every energy pair (conf1, conf2) one equation is formed:

$$E_2 - E_1 = \text{factor} \sum_{\langle ij \rangle} \sum_{\alpha\beta} J_{ij}^{\alpha\beta}
\left( S_{i\alpha}^{(2)} S_{j\beta}^{(2)} - S_{i\alpha}^{(1)} S_{j\beta}^{(1)} \right)$$

where the spin components of each configuration are rotated into the frame of pair (*i*, *j*).
With *n* pairs there are 9*n* unknowns, ordered `pair-1_xx, pair-1_xy, …, pair-n_zz`.

---

## Directory structure

```
anisotropic_exchange/
├── main.jl                             # fits the J matrices from input.yaml
├── input.yaml                          # example input (3 spin-3/2 sites, 3 pairs, 27 energy pairs)
├── get_emat.py                         # builds a pair's local reference frame from a structure
├── format_energies.py                  # converts energies.dat rows into input.yaml energy pairs
├── get_average_energies.py             # averages the two halves of energies.dat
└── transform_J_matrix_local2global.py  # rotates a fitted J matrix into the global frame
```

---

## `input.yaml`

```yaml
spins:                 # spin quantum number of each site
  - 1.5
  - 1.5
  - 1.5

exchange:              # one entry per pair; defines the unknown J matrices
  - pair: [1, 2]
    reference_frame: [ex_x, ex_y, ex_z,  ey_x, ey_y, ey_z,  ez_x, ez_y, ez_z]
  ...

energies:              # one entry per equation
  - pair: [1, 2]       # frame in which the angles below are given
    conf1: [theta_1, phi_1, theta_2, phi_2, ..., energy]
    conf2: [theta_1, phi_1, theta_2, phi_2, ..., energy]
  ...
```

- `reference_frame` lists the local axes **e**<sub>x</sub>, **e**<sub>y</sub>, **e**<sub>z</sub> (in
  global coordinates) one after another. This is the format printed by `get_emat.py`.
- In `energies`, each configuration gives polar and azimuthal angles (degrees) for every site,
  followed by the total energy in eV. The angles are given in the reference frame of the
  listed `pair`.
- At least 9 energy pairs are needed per exchange pair. With exactly that many, the system is
  solved exactly. With more, it is solved by least squares and the residual is printed.

---

## `main.jl`

| Function | Description |
|----------|-------------|
| `sph2cart(r, theta, phi, deg)` | Spherical → Cartesian spin vector |
| `transform_reference_frame(...)` | Rotates spin vectors from one pair frame to another |
| `fatten_one_energy_pair(...)` | Expresses both configurations of one energy pair in every pair frame |
| `get_one_row_of_coefficients(...)` | Builds one row of the coefficient matrix |
| `get_coefficient_matrix(..., factor)` | Builds the full coefficient matrix and the energy differences (meV) |
| `get_Js(..., factor, quantum_correction)` | Solves for the J components (exact or least-squares; errors if underdetermined); optionally applies the empirical correction 2*S*<sub>*i*</sub>*S*<sub>*j*</sub> / (2*S*<sub>*i*</sub>*S*<sub>*j*</sub> + min(*S*<sub>*i*</sub>, *S*<sub>*j*</sub>)) (Ruiz 2003) |
| `print_Js(n_spin_pair, Js, iformat)` | `iformat = 1`: labelled per pair; `iformat = 2`: one line of all components |
| `main(dirname)` | Reads `dirname/input.yaml`, fits, and prints the J components |

### Usage

```bash
julia main.jl                 # uses input.yaml in the current directory
julia main.jl path/to/dir     # uses path/to/dir/input.yaml
```

`factor` and `quantum_correction` are set in `main`.

---

## Helper scripts

| Script | Description |
|--------|-------------|
| `get_emat.py` | Reads `left_hand_relaxed.vasp` and builds the frame of pair (*i*, *j*): **e**<sub>x</sub> along the *i* → *j* bond, **e**<sub>z</sub> = the part of the global *z* axis perpendicular to the bond, **e**<sub>y</sub> = **e**<sub>z</sub> × **e**<sub>x</sub>. Prints it as a `reference_frame:` line or as a 3×3 array |
| `format_energies.py` | Reads `energies.dat` (columns: θ, φ per site, energy) and prints `input.yaml` energy pairs for one site pair, taking rows *k* and *k* + 9 as conf1 and conf2 |
| `get_average_energies.py` | Averages the energies of row *k* and row *k* + *m*/2 of `energies.dat`; writes `energies_average.dat` |
| `transform_J_matrix_local2global.py` | Rotates a J matrix from a local frame to the global frame, here for the three C₃-related pairs (0°, 120°, 240° about *z*) |

---

## Dependencies

- Julia with `YAML.jl` (`LinearAlgebra` and `Printf` are standard libraries)
- Python ≥ 3.8 with `numpy`; `ase` for `get_emat.py`; `scipy` for `transform_J_matrix_local2global.py`
