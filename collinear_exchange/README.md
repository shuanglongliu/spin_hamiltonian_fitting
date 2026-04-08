# collinear_exchange

Extract isotropic exchange coupling constants *J* from collinear broken-symmetry DFT calculations.

The workflow has three stages, each handled by a dedicated directory:

1. **`orca/`** – run broken-symmetry ORCA calculations on a molecular cluster
2. **`vasp/`** – run broken-symmetry VASP calculations on a periodic or cluster system
3. **`fit_J/`** – fit the *J* values from the collected DFT energies

---

## Spin Hamiltonian

$$H = \text{factor} \times \sum_k J_k \sum_{\langle ij \rangle \in k} \mathbf{S}_i \cdot \mathbf{S}_j$$

`factor` is typically `−2` (Heisenberg convention). The *J* values are reported in cm⁻¹.

---

## Directory structure

```
collinear_exchange/
├── orca/
│   ├── data.py       # spin configs, ORCA input template, SLURM script
│   ├── main.py       # BsDft / BsDftManager classes
├── vasp/
│   ├── data.py       # INCAR, POSCAR, KPOINTS templates, SLURM scripts
│   └── main.py       # helper functions for VASP workflow
└── fit_J/
    ├── get_Js.py     # reads energies.dat + pairs.yaml, solves for Js
    ├── energies.dat  # output of orca/main.py or vasp/main.py (auto-generated)
    └── pairs.yaml    # user-provided spin-pair topology
```

---

## `orca/` — ORCA broken-symmetry calculations

Target system in the current example: **Nd₃Mn₈** (11 magnetic centers).

### `data.py`

| Symbol | Description |
|--------|-------------|
| `spin_only_Szs` | Per-site spin-only \|*m*s\| values (length = number of magnetic sites) |
| `effective_spins` | Per-site effective spin quantum numbers used when writing `energies.dat` |
| `expected_magmoms` | Reference Löwdin spin populations from the high-spin calculation, used for sanity checks |
| `configs` | List of spin configuration strings (`'u'` = spin-up, `'d'` = spin-down); includes all-up HS reference, single-site and two-site spin flips |
| `input_file` | ORCA input template (TPSSh / ZORA-TZVP; `{flip_spin}`, `{final_ms}`, `{multiplicity}` placeholders) |
| `job_script` | SLURM batch script for the HPC cluster |

### `main.py`

#### `BsDft` — single configuration

```python
BsDft(config, spin_only_Szs, restart=False)
```

| Method | Description |
|--------|-------------|
| `create_work_dir()` | Creates a subdirectory named after the config string |
| `get_input_file()` | Writes `orca.inp`, `orca.job`; symlinks `orca.xyz` and the reference `dft.gbw` |
| `submit_job()` | `sbatch orca.job` from the working directory |
| `check_convergence()` | Counts `"SCF CONVERGED"` in `orca.out` (expects 2 for initial spin-flip run, 1 for restart) |
| `backup_files()` | Removes scratch files; archives `orca.inp/out/gbw` → `run_N.*`; updates `dft.gbw` symlink |
| `resubmit_job()` | Resubmits if not converged |
| `get_magnetic_moments()` | Parses Löwdin spin populations from `orca.out`; validates sign and magnitude; writes `magnetic_moments.csv` |
| `get_energy()` | Returns the final SCF energy (Ha) after convergence and magnetic-moment checks |

`restart=True` comments out `FlipSpin`/`FinalMs` directives so ORCA reads the spin state from the existing GBW file.

#### `BsDftManager` — all configurations

```python
BsDftManager(configs, spin_only_Szs, restart=False)
```

Iterates `BsDft` over a list of configurations and provides the same methods in bulk. `get_energies()` collects all converged energies and writes `energies.dat`.

#### Typical workflow

```python
# orca/main.py __main__ block – uncomment steps as needed
manager = BsDftManager(configs=configs, spin_only_Szs=spin_only_Szs, restart=False)
manager.create_work_dirs()    # Step 1
manager.get_input_files()     # Step 2
manager.submit_jobs()         # Step 3
manager.check_convergence()   # Step 4a
manager.backup_files()        # Step 4b – backs up non-converged runs
manager.resubmit_jobs()       # Step 5
manager.get_energies()        # Step 6 – writes energies.dat
```

---

## `vasp/` — VASP broken-symmetry calculations

Target system in the current example: **Ni₇** (7 magnetic centers in a molecular cluster).

### `data.py`

| Symbol | Description |
|--------|-------------|
| `incar` | INCAR template (HSE06 functional; `{magmom}` and `{nupdown}` placeholders; comments for R2SCAN, B3LYP, LDA+U, vdW alternatives) |
| `poscar` | POSCAR for the Ni₇O₂₆C₈₂H₁₆₆ cluster in a large box |
| `kpoints` | Gamma-point KPOINTS |
| `job_script_cpu` | SLURM script for CPU nodes (NERSC Perlmutter, `vasp_gam`) |
| `job_script_gpu` | SLURM script for GPU nodes with automatic STOPCAR after 22 h |
| `job_script` | Alias pointing to the active script (currently `job_script_gpu`) |

### `main.py`

| Function | Description |
|----------|-------------|
| `create_dir(directories)` | Creates one subdirectory per configuration string |
| `dir2magmom(dirname, magmoms)` | Converts a config string + per-site moments to a VASP `MAGMOM` string (`+` for `'u'`, `−` for `'d'`) |
| `get_total_magmom(dirname, magmoms)` | Returns the signed total magnetic moment (`NUPDOWN`) |
| `get_inp(directories, magmoms)` | Writes `INCAR`, `POSCAR`, `KPOINTS`, `vasp.job`; symlinks `POTCAR` |
| `submit_job(directories)` | `sbatch vasp.job` for each directory |
| `check_convergences(directories)` | Greps `OUTCAR` for `"EDIFF is reached"` |
| `check_magmoms(directories)` | Reads site magnetic moments from the last `"magnetization (x)"` block in `OUTCAR` |
| `get_energies(directories)` | Reads total energies (eV, `"sigma"` line in `OUTCAR`); writes `energies.dat` |

#### Typical workflow

```python
# vasp/main.py __main__ block – uncomment steps as needed
directories = ["uuuuuuu", "uduuuuu", ...]
magmoms     = [2, 2, 2, 2, 2, 2, 2]
effective_spins = [1, 1, 1, 1, 1, 1, 1]

create_dir(directories)
get_inp(directories, magmoms)
submit_job(directories)
check_convergences(directories)
check_magmoms(directories)
get_energies(directories)        # writes energies.dat
```

---

## `fit_J/` — Extract exchange coupling constants

### Input files

**`energies.dat`** (generated by `orca/main.py` or `vasp/main.py`)
```
# S1    S2  ...  Energy (eV or Ha)
  1.0   1.0 ...  -123.456789
 -1.0   1.0 ...  -123.451234
  ...
```
The first header line must contain the energy unit (`eV` or `Ha`).

**`pairs.yaml`** (user-provided)
```yaml
J1: [[1, 2], [3, 4]]   # 1-indexed site pairs sharing the same J
J2: [[1, 3]]
...
```

### `get_Js.py`

| Function | Description |
|----------|-------------|
| `read_data(energy_threshold)` | Loads `energies.dat`; converts energies to cm⁻¹ relative to the first configuration; optionally removes configurations more than `energy_threshold` cm⁻¹ above the minimum |
| `read_pairs()` | Loads `pairs.yaml`; returns a dict mapping J labels to lists of 1-indexed site-pair lists |
| `get_sumSS_vector(data, pairs)` | Computes Σ *Si*·*Sj* for one *J* over all configurations |
| `get_sumSS_matrix(data, pairs_dict)` | Builds the (n_J × n_config) coefficient matrix |
| `get_Js(data, sumSS_matrix, factor)` | Solves `factor × C × J = ΔE`; exact solve for square systems, least-squares for overdetermined; writes `Js.csv`; saves `fit_quality.png` when overdetermined |

#### Usage

```python
# fit_J/get_Js.py __main__ block
data        = read_data(energy_threshold=1000)   # exclude outliers > 1000 cm^-1
pairs_dict  = read_pairs()
sumSS_matrix = get_sumSS_matrix(data, pairs_dict)
Js          = get_Js(data, sumSS_matrix, factor=-2)
```

Results are saved to `Js.csv`:
```
# Js in cm^-1, factor = -2
J1,J2,...
12.34,-5.67,...
```

---

## Dependencies

- Python ≥ 3.8
- `numpy`
- `pyyaml`
- `matplotlib` (only for the overdetermined least-squares fit plot)
- ORCA ≥ 6 (for `orca/`) and/or VASP ≥ 6.4 (for `vasp/`)
- SLURM (for HPC job submission)
