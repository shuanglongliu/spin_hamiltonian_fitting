# Magnetic Anisotropy

Python codes for managing VASP calculations of magnetic anisotropy. The codes in the global directory are for collinear spin alignments. The codes in the local directory are for arbitrary spin alignments.

Usage: 

Step 1: Modify the vairables in data.py

Step 2: Modify main.py to set up VASP calculations.

Step 3: Modify main.py to submit jobs using Slurm.

Step 4: Modify main.py to collect results.


Requirement: A POTCAR file in the root directory where you set up all DFT calculations.

The spherical-coordinate and reference-frame helpers used by the scripts in `data_processing/`
are in `data_processing/geometry.py`.
