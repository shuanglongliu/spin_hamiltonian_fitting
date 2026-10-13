# Magnetic Anisotropy

Python codes for managing VASP calculations of magnetic anisotropy. The codes in the global directory are for collinear spin alignments. The codes in the local directory are for arbitrary spin alignments.

Usage: 

Step 1: Modify the vairables in data.py

Step 2: Modify main.py to set up VASP calculations.

Step 3: Modify main.py to submit jobs using Slurm.

Step 4: Modify main.py to collect results.


Requirement: A POTCAR file in the root directory where you set up all DFT calculations.

The fitting scripts `fit_D_and_E.py`, `fit_D_tensor.py`, `fit_energy_surface.py` and
`get_average_energies.py` in `data_processing/` import `pylib` from
[shuanglongliu/programs](https://github.com/shuanglongliu/programs) (`python/pylib.py`); add that
directory to `PYTHONPATH`.
