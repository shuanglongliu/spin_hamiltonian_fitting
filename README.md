# spin_hamiltonian_fitting

Tools for obtaining spin Hamiltonian parameters from DFT calculations and fitting.

| Directory | Parameters | Description |
|---|---|---|
| [`collinear_exchange/`](collinear_exchange/) | Isotropic exchange *J* | Broken-symmetry ORCA/VASP calculations and fitting of *J* from collinear DFT energies |
| [`magnetic_anisotropy/`](magnetic_anisotropy/) | Single-ion anisotropy (*D*, *E*, **D** tensor) | VASP workflows for global (collinear) and local (arbitrary) spin orientations, plus fitting scripts |
| [`polarization/`](polarization/) | Spin-driven electric polarization | Fitting the coefficients of the dot- and cross-product spin terms in the polarization |

Each subdirectory was previously its own repository; its full commit history is preserved here
(`git log -- <directory>/`).

## License

Each subdirectory keeps its original license:

- `collinear_exchange/`, `polarization/`: GNU AGPL v3
- `magnetic_anisotropy/`: GNU GPL v3
