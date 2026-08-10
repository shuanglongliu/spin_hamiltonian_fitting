import os

root_dir = os.path.dirname(os.path.abspath(__file__))  + "/"

incar = """
SYSTEM = vasp

#### symmetry ####
#ISYM = 0

#### system size ####
#NBANDS = 448
#NELECT = 664.0

#### accuracy ####
PREC = Accurate
ENCUT = 500
LREAL = A

#### parallelization ####
#KPAR = 4
NCORE = 8

#### density functional ####

#=== LDA/GGA ===#
#GGA = PE

#=== METAGGA ===#
#METAGGA = R2SCAN

#=== Hybrid: HSE06 ===#
LHFCALC = T 
GGA = PE
HFSCREEN = 0.2 
PRECFOCK = Accurate

#=== vdW-DF: optB86b ===#
#GGA      = MK 
#PARAM1   = 0.1234 
#PARAM2   = 1.0
#AGGAC    = 0.0
#LUSE_VDW = .TRUE.
#LASPH    = .TRUE.

#=== LIBXC ===#
#GGA = LIBXC
#LIBXC1 = GGA_X_PBE
#LIBXC2 = GGA_C_PBE

#### empirical vdW ####
#IVDW = 11

#### LDA+U ####
#LDAU = T
#LDAUTYPE = 2
#LDAUPRINT = 1
#LDAUL =   2   -1    -1   -1   -1   -1   -1  
#LDAUU =   2.5  0.0   0.0  0.0  0.0  0.0  0.0
#LDAUJ =   0.0  0.0   0.0  0.0  0.0  0.0  0.0

#### electronic optimization ####
EDIFF = 1E-6
NELM = 120
#NELMIN = 10
NELMDL = -12
AMIX = 0.3
AMIX_MAG = 0.6
AMIX_MIN = 0.1
BMIX = 0.0001
BMIX_MAG = 0.0001
ALGO = All
# TIME = 0.35

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
SAXIS = 0 0 1
MAGMOM = {mstring:s} 3000*0.0

#### magnetism: constraint ####
I_CONSTRAINED_M = 4
M_CONSTR = {mstring:s} 3000*0.0

LAMBDA =  15.0
RWIGS = {rwigs:s}

#### magnetism: orbital moment ####
LORBMOM = T

#### charge, wavefunction ####
ISTART = 0
ICHARG = 1
LWAVE = T
LCHARG = F
LAECHG = F
LMAXMIX = 6

#### dos ####
ISMEAR = 0
SIGMA = 0.01
#NEDOS = 2501
#EMIN = -15
#EMAX = 10
LORBIT = 10

#### Partial Charge ####
#LPARD = T
#EINT = -5.0000  -4.7355
#EINT = -4.7355 -4.0500

#### wann ####
#LWANNIER90 = .T.
#LWRITE_UNK = .TRUE.

#### polarization ####
#IDIPOL = 4
#LDIPOL = T
#DIPOL = 0.5 0.5 0.5
#LCALCPOL = T

#### occupation matrix control ####
#OCCEXT = 1
"""

kpoints = """kpoints
0
gamma
 1   1   1
 0.0 0.0 0.0
"""

poscar = """
"""

job_script = """#!/bin/bash

#SBATCH -A m3346
#SBATCH -J Mn7
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH -t 12:00:00
#SBATCH -N 2
#SBATCH -n 256
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:1704659

module load vasp/6.3.2-cpu

export OMP_NUM_THREADS=1

source ~/.bash_aliases; keep_log; # vp_continue_relaxation

srun -n 256 -c 2 --cpu-bind=cores vasp_ncl
"""

job_script = """#!/bin/bash

#SBATCH --account=m2qm-efrc
#SBATCH --qos=m2qm-efrc-b
#SBATCH --job-name=Gd2CeMn8
#SBATCH --partition=hpg-default
#SBATCH --constraint=el9
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=512gb
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mail-type=None
#SBATCH --mail-user=shuan.liu@northeastern.edu
#SBATCH --error=error
#SBATCH --output=output

module purge
module load intel/2025.1.0 openmpi/5.0.7 vasp/6.4.3

srun --mpi=pmix_v5 vasp_ncl
"""

emats_file = []

for i in range(100):
    emats_file.append( [[1,0,0],[0,1,0],[0,0,1]] )

