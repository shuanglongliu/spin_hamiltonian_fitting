# spin_hamiltonian_fitting

Spin Hamiltonian parameters from DFT: set up and harvest VASP/ORCA jobs for chosen spin
configurations, then fit the parameters to the energies or dipoles.

| Parameters | DFT jobs (template) | Fit | Docs |
|---|---|---|---|
| Isotropic exchange *J* | `collinear_vasp`, `collinear_orca` | `shf isotropic-exchange` | [isotropic_exchange.md](docs/isotropic_exchange.md) |
| Exchange matrices **J** (3×3 per pair) | `site_vasp` | `shf anisotropic-exchange fit` | [anisotropic_exchange.md](docs/anisotropic_exchange.md) |
| Single-ion anisotropy (*D*, *E*, **D** tensor) | `saxis_vasp`, `site_vasp` | `shf anisotropy <task>` | [single_ion_anisotropy.md](docs/single_ion_anisotropy.md) |
| Spin-driven electric polarization | `site_vasp` (`dipoles`) | `shf polarization <model>` | [polarization.md](docs/polarization.md) |

The DFT job templates are described in [dft_jobs.md](docs/dft_jobs.md).

## Install

```bash
pip install -e .            # Python >= 3.8; numpy, scipy, pyyaml, matplotlib
pip install -e ".[structure]"   # also ase, for `shf anisotropic-exchange frame`
```

This installs the `shf` package and command. Job scripts (`main.py`) import `shf`, so they run
from any directory.

## A new system

```bash
shf new site_vasp ~/runs/Mn3      # copy a template: data.py + main.py
cd ~/runs/Mn3
# 1. data.py: INCAR, POSCAR, KPOINTS, job script; put POTCAR here
# 2. main.py: the spin configurations
python main.py print              # check the configurations and directory names
python main.py setup              # write the inputs
python main.py submit
python main.py check              # later: convergence
python main.py energies           # energies.dat (dipoles: dipoles.dat)
```

then fit, e.g. `shf new fit_anisotropy ~/runs/Mn3/fit`, edit `anisotropy.yaml`, and run
`shf anisotropy fit-DE`. `shf -h`, `shf <command> -h` and `python main.py help` list the options.

## Layout

```
src/shf/
├── geometry.py, sphere.py        # spherical coordinates, frames, Euler rotations; direction grid
├── spin_operators.py, units.py
├── linfit.py, tables.py, plotting.py
├── dft/                          # job sets: VASP (collinear, SAXIS, per-site MAGMOM), ORCA, Slurm
├── models/                       # isotropic and anisotropic exchange, single-ion anisotropy, polarization
└── cli.py                        # the shf command
templates/                        # starting points for `shf new`
├── collinear_vasp, collinear_orca, saxis_vasp, site_vasp      # data.py + main.py
└── fit_*                         # fit inputs with example data
docs/
tests/                            # python -m unittest discover tests
```

Each subdirectory of the earlier layout (`collinear_exchange/`, `anisotropic_exchange/`,
`magnetic_anisotropy/`, `polarization/`) was its own repository; the full history is in
`git log`.

## License

GNU General Public License v3.0; see [LICENSE](LICENSE).
