# spin_hamiltonian_fitting

Tools for obtaining spin Hamiltonian parameters from DFT calculations and fitting.

| Directory | Parameters | Description |
|---|---|---|
| [`collinear_exchange/`](collinear_exchange/) | Isotropic exchange *J* | Broken-symmetry ORCA/VASP calculations and fitting of *J* from collinear DFT energies |
| [`anisotropic_exchange/`](anisotropic_exchange/) | Exchange matrices **J** (3×3 per pair) | Fitting full exchange coupling matrices from noncollinear DFT energies in pair-local reference frames |
| [`magnetic_anisotropy/`](magnetic_anisotropy/) | Single-ion anisotropy (*D*, *E*, **D** tensor) | VASP workflows for global (collinear) and local (arbitrary) spin orientations, plus fitting scripts |
| [`polarization/`](polarization/) | Spin-driven electric polarization | Fitting the coefficients of the dot- and cross-product spin terms in the polarization |

Each subdirectory was previously its own repository; its full commit history is preserved here
(`git log -- <directory>/`).

## License

GNU General Public License v3.0; see [LICENSE](LICENSE).
