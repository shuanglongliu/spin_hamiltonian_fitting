# Magnetic Anisotropy

Python codes for managing VASP calculations of magnetic anisotropy. The codes in the global directory are for collinear spin alignments. The codes in the local directory are for arbitrary spin alignments.

Usage: 

Step 1: Modify the vairables in data.py

Step 2: Modify main.py to set up VASP calculations.

Step 3: Modify main.py to submit jobs using Slurm.

Step 4: Modify main.py to collect results.


Requirement: A POTCAR file in the root directory where you set up all DFT calculations.

`base.py` holds the helpers shared by every directory here: spherical/cartesian conversion,
reference-frame changes, the discretized `sphere`, and the `magmoms` class used to write MAGMOM.
The job directories (`global_magnetic_anisotropy/`, `local_magnetic_anisotropy/*/`) and the
scripts in `data_processing/` import it from `magnetic_anisotropy/`, so keep this directory tree
together (run in the repository checkout, or copy the whole `magnetic_anisotropy/` folder).
