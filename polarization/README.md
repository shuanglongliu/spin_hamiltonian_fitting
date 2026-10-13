# polarization

Fit the spin-driven electric polarization of a magnetic cluster to DFT dipoles.

Each fit reads `input.yaml` and `dipoles.dat` from the working directory and minimizes
|*A* *x* − *y*|², where *y* holds the change in the dipole of every spin configuration relative to
the first one. When the data do not determine every parameter, the scripts return the
minimum-norm solution and print a warning that says how many parameter combinations are
undetermined.

```
polarization/
├── alpha_and_pi/          # scalar alpha (cross-product term) and vector pi (dot-product term) per pair
│   ├── main.py
│   ├── input.yaml
│   └── dipoles.dat
└── full_p_matrix/         # full 3×3×3 polarization tensor per pair
    ├── notes.txt
    ├── version_1/         # all pairs as unknowns, symmetry imposed as constraints
    └── version_2/         # only irreducible pairs as unknowns, others obtained by rotation
```

Usage, in any of the three directories:

```bash
python main.py                 # uses input.yaml and dipoles.dat in the current directory
python main.py path/to/dir     # uses path/to/dir/
```

---

## `dipoles.dat`

One row per spin configuration: θ, φ (degrees) for every site, then *p*<sub>x</sub>, *p*<sub>y</sub>,
*p*<sub>z</sub> (e Å). The first row is the reference configuration.

---

## `alpha_and_pi/`

$$\mathbf{p} = \sum_{\langle ij \rangle} \Big[ \alpha_{ij} \big( \mathbf{S}_i (\mathbf{S}_j \cdot \mathbf{e}_{ij}) - (\mathbf{S}_i \cdot \mathbf{e}_{ij}) \mathbf{S}_j \big) + \boldsymbol{\pi}_{ij} (\mathbf{S}_i \cdot \mathbf{S}_j) \Big]$$

**e**<sub>*ij*</sub> is the unit vector from site *i* to site *j*.

`input.yaml`:

| Key | Description |
|-----|-------------|
| `spins` | Spin quantum number of each site |
| `pairs` | 1-indexed site pairs |
| `positions` | Cartesian site positions (Å), used for **e**<sub>*ij*</sub> |
| `constraints_alpha` | `[pair_i, pair_j]`: α(pair_i) = α(pair_j) |
| `constraints_pi` | `[pair_i, pair_j, α, β, γ]`: **π**(pair_j) = *R* **π**(pair_i), with *R* from extrinsic zyz Euler angles (degrees) |

Outputs: `log` (α and **π** per pair, maximal and average deviation) and `fitting_quality.pdf`.

---

## `full_p_matrix/`

$$p_\gamma = \sum_{\langle ij \rangle} \sum_{\alpha\beta} S_{i\alpha}\, P^{ij}_{\gamma\alpha\beta}\, S_{j\beta}$$

Elements are stored as `x[27 * k + 9 * gamma + 3 * alpha + beta]` (β fastest, then α, γ, and
the pair *k*).

### `version_1/`

All pairs are unknowns. Each `constraints` entry `[pair_i, pair_j, α, β, γ]` imposes
*P*<sup>pair_j</sup> = (*R* ⊗ *R* ⊗ *R*) *P*<sup>pair_i</sup>.
Outputs: `P_total.log`, `fit.log` and `fitting_quality.pdf`.

### `version_2/`

Only the pairs listed in `irr_pairs` are unknowns. Each `transformations` entry
`[pair_i, pair_j, α, β, γ]` builds *P*<sup>pair_j</sup> from the irreducible pair_i.
`emats` gives the local reference frame of each irreducible pair (**e**<sub>x</sub>,
**e**<sub>y</sub>, **e**<sub>z</sub> one after another).

By default `main` writes `P_for_latex.log`: the tensors in the local frames, split into
isotropic (*P*<sup>J</sup>), antisymmetric (*P*<sup>D</sup>) and anisotropic symmetric
(*P*<sup>Γ</sup>) parts, in units of 10⁻⁵ e Å. Uncomment the calls in `main` for:

| Function | Output |
|----------|--------|
| `check_fitting_quality` | `fit.log` (absolute and percentage deviations), `fitting_quality.pdf` |
| `save_results(x, inp, modifier)` | `P_total`, `P_J`, `P_D`, `P_Gamma` and `pi_and_d` logs (global frame with `""`, local frames with `"_loc"`) |
| `get_relative_contribution` | `relative_contribution.log`: share of each irreducible pair and of its J, D, Γ parts |

With the same data, version 1 and version 2 give the same result; see `full_p_matrix/notes.txt`.

---

## Dependencies

- Python ≥ 3.8 with `numpy`, `scipy`, `pyyaml` and `matplotlib`
