"""Static templates and paths used to drive VASP runs.

Edit ``incar``, ``poscar``, ``kpoints``, and ``job_script`` here to change
the calculation; ``common.py`` writes them verbatim into each run directory.
"""

import os


# Directory layout
root_dir = os.path.dirname(os.path.abspath(__file__))
wave_dir = "0"  # Reference run that supplies WAVECAR to subsequent jobs (the other branch used "6476401").


incar = """
SYSTEM = vasp

#### sym ####
#ISYM = 0

#### system size ####
#NBANDS = 448
#NELECT = 383.0

#### accuracy ####
PREC = Accurate
ENCUT = 450
LREAL = F

#### parallelization ####
# KPAR = 4
NCORE = 8

#### electronic optimization ####
EDIFF = 1E-8
NELM = 90
#NELMIN = 10
#NELMDL = -8
AMIX = 0.2
AMIX_MAG = 0.4
AMIX_MIN = 0.05
BMIX = 0.0001
BMIX_MAG = 0.0001
ALGO = All

#### structural relaxation ####
#NSW = 500
#IBRION = 2
#ISIF = 2
#EDIFFG = -0.02
#POTIM = 0.3

#### magnetism: accuracy ####
LASPH = T
GGA_COMPAT = F

#### magnetism: collinear spin ####
#ISPIN = 2
#MAGMOM = 3.0 3.0 3.0 1000*0.0
#NUPDOWN = 9

#### magnetism: noncollinear spin, SOC ####
LSORBIT = T
SAXIS = {sstring:s}
MAGMOM = 0.0 0.0 5.0 \\
         300*0.0

#### magnetism: constraint ####
I_CONSTRAINED_M = 1
M_CONSTR =  0.0 0.0 5.0 \\
         300*0.0

LAMBDA =  10.0
RWIGS = 1.323 1.524 0.820

#### magnetism: orbital moment ####
LORBMOM = T

#### charge, wavefunction ####
ISTART = 1
ICHARG = 0
LWAVE = F
LCHARG = F
LAECHG = F
LMAXMIX = 4

#### dos ####
ISMEAR = 0
SIGMA = 0.1
# NEDOS = 501
# EMIN = -15
# EMAX = 10
# LORBIT = 11

#### Partial Charge ####
#LPARD = T
#EINT = -5.0000  -4.7355
#EINT = -4.7355 -4.0500

#### vdW ####
# IVDW = 11

#### LDA+U ####
# LDAU = T
# LDAUTYPE = 2
# LDAUPRINT = 1
# LDAUL =   2    2   -1
# LDAUU =   5.5  3.0  0.0
# LDAUJ =   0.0  0.0  0.0

#### HSE ####
#LHFCALC = T
#HFSCREEN = 0.2
#PRECFOCK = Accurate
#ALGO = All
#TIME = 0.35

#### wann ####
#LWANNIER90 = .T.
#LWRITE_UNK = .TRUE.

### polarization ###
#EFIELD_PEAD = 0.000 0.000 0.000
#IDIPOL = 3
#LMONO = T
#LDIPOL = T
#LCALCPOL = T
#DIPOL = 0.5 0.5 0.5

### occupation matrix control ###
#OCCEXT = 1
"""

kpoints = """kpoints
0
gamma
 15  15  15
 0.0 0.0 0.0
"""

poscar = """Mn Mg3 O4 - experiment lattice constant
1.0
   4.2112    0.0000    0.0000
   0.0000    4.2112    0.0000
   0.0000    0.0000    4.2112
Mn Mg O
1 3 4
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Mg2+
   0.0000000000000000    0.5000000000000000    0.5000000000000000 Mg2+
   0.5000000000000000    0.0000000000000000    0.5000000000000000 Mg2+
   0.5000000000000000    0.5000000000000000    0.0000000000000000 Mg2+
   0.0000000000000000    0.0000000000000000    0.5000000000000000 O2-
   0.0000000000000000    0.5000000000000000    0.0000000000000000 O2-
   0.5000000000000000    0.0000000000000000    0.0000000000000000 O2-
   0.5000000000000000    0.5000000000000000    0.5000000000000000 O2-
"""

# Slurm submission script for Perlmutter CPU partition.
# (A GPU variant is kept below as a reference; switch by reassigning
# ``job_script`` to the desired template.)
job_script_gpu = """#!/bin/bash
#SBATCH -A m3346_g
#SBATCH -J Cr3
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -G 4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --ntasks-per-node=4
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45514509

module load vasp/6.4.3-gpu

source ~/.bash_aliases; keep_log

srun -n 4 -c 32 --cpu-bind=cores -G 4 --gpu-bind=none vasp_ncl
"""

job_script = """#!/bin/bash

#SBATCH -A m3346
#SBATCH -J Mn-MgO
#SBATCH -C cpu
#SBATCH --qos=shared
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 2
#SBATCH --mem=128GB
#SBATCH --output=output
#SBATCH --error=error
#SBATCH --mail-type=None
#SBATCH --mail-user=shuan.liu@northeastern.edu

source ~/.bash_aliases; keep_log

module load vasp/6.4.3-cpu
export OMP_NUM_THREADS=1

srun -n 64 -c 2 --cpu-bind=cores vasp_ncl
"""


# Default identity transformation: the local frame coincides with the global
# frame. ``vasp_job_ncl.setup_local_ref_frames(from_file=True)`` reads this.
emat_file = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
]
