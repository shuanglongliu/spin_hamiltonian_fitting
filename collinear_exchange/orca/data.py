"""
Static data for broken-symmetry DFT calculations on the Nd3Mn8 cluster.

Defines spin configurations, per-site spin quantum numbers, ORCA input
file templates, and the SLURM job script used by main.py.
"""

import os

# Paths anchored to this file's location
root_dir = os.path.dirname(os.path.abspath(__file__))

# High-spin (all-up) reference calculation used as MO starting guess
hs_dir = os.path.join(root_dir, "uuuu")

# Per-site spin quantum numbers |ms| for the 11 magnetic centers:
spin_only_Szs = [0.5, 0.5, 0.5, 0.5]

# Effective spin quantum numbers
effective_spins = [0.5, 0.5, 2.5, 2.5]

# Expected local magnetic moments (in μB) for the HS reference, from Loewdin analysis:
expected_magmoms = [
1.0,
1.0,
1.0,
1.0,
]

# Broken-symmetry spin configurations 
#   'u' = spin-up (ms = +|ms|),  'd' = spin-down (ms = −|ms|)
configs = [
    "uuuu",
    "uddd",
    "uddu",
    "udud",
    "uduu",
    "uudd",
    "uudu",
    "uuud",
]

configs = [ "uuuu" ]

n_configs = len(configs)

n_cpus = 32
mem_per_core = 16000

# ─────────────────────────────────────────────────────────────────────────────
# ORCA input file templates
# ─────────────────────────────────────────────────────────────────────────────

input_file_hs = """
# --- Method & relativistic treatment ---
# ZORA scalar relativistics are essential for lanthanides
# TPSSh gives reliable J values for transition/lanthanide systems
! PBE ZORA DefGrid2 TightSCF MORead KDIIS NoTRAH # SlowConv

%basis
    # All-electron ZORA basis for lanthanide atoms — must use SARC
    NewGTO Ce "SARC-ZORA-TZVP" end

    # ZORA-recontracted triple-zeta for transition metal and light atoms
    NewGTO Cu "ZORA-def2-TZVP"  end
    NewGTO O  "ZORA-def2-SVP"  end
    NewGTO N  "ZORA-def2-SVP"  end
    NewGTO C  "ZORA-def2-SVP"  end
    NewGTO H  "ZORA-def2-SVP"   end   # SVP for H is sufficient
end

# --- Parallelization ---
%pal nprocs {n_cpus:d} end

# --- Memory per core in MB ---
%maxcore {mem_per_core:d}

# --- SCF options ---
%scf
    MaxIter   500

    # Use larger DIIS to help convergence with open-shell lanthanide
    DIISMaxEq 15

    # Trusted region augmentation Hessian (TRAH): 2nd order SCF
    # AutoTRAH  true

    # Use the converged HS orbitals as starting guess
    MOInp "dft.gbw"

    # Flip spins on magnetic centers (0-based atom indices)
    # {flip_spin_tag:s} {flip_spin:s}
    # {final_ms_tag:s} {final_ms:.1f}
end

# --- Output: natural population analysis useful for spin analysis ---
%output
    Print[P_Mulliken]   1
    Print[P_Loewdin]  1
end

# ================================================================
#  CHARGE and MULTIPLICITY
# ================================================================
* xyzfile  0  {multiplicity:d}   orca.xyz
"""

input_file_bs = """
# --- Method & relativistic treatment ---
# ZORA scalar relativistics are essential for lanthanides
# TPSSh gives reliable J values for transition/lanthanide systems
! TPSSh ZORA DefGrid3 TightSCF KDIIS NoTRAH # SlowConv MORead

%basis
    # All-electron ZORA basis for lanthanide atoms — must use SARC
    NewGTO Ce "SARC-ZORA-TZVP" end

    # ZORA-recontracted triple-zeta for transition metal and light atoms
    NewGTO Cu "ZORA-def2-TZVP"  end
    NewGTO O  "ZORA-def2-TZVP"  end
    NewGTO N  "ZORA-def2-TZVP"  end
    NewGTO C  "ZORA-def2-TZVP"  end
    NewGTO H  "ZORA-def2-SVP"   end   # SVP for H is sufficient
end

# --- Parallelization ---
%pal nprocs {n_cpus:d} end

# --- Memory per core in MB ---
%maxcore {mem_per_core:d}

# --- SCF options ---
%scf
    MaxIter   500

    # Use larger DIIS to help convergence with open-shell lanthanide
    DIISMaxEq 15

    # Trusted region augmentation Hessian (TRAH): 2nd order SCF
    # AutoTRAH  true

    # Use the converged HS orbitals as starting guess
    MOInp "dft.gbw"

    # Flip spins on magnetic centers (0-based atom indices)
    {flip_spin_tag:s} {flip_spin:s}
    {final_ms_tag:s} {final_ms:.1f}
end

# --- Output: natural population analysis useful for spin analysis ---
%output
    Print[P_Mulliken]   1
    Print[P_Loewdin]  1
end

# ================================================================
#  CHARGE and MULTIPLICITY
# ================================================================
* xyzfile  0  {multiplicity:d}   orca.xyz
"""

job_script = """#!/bin/bash

#SBATCH -A m2qm-efrc
#SBATCH --qos=m2qm-efrc-b
#SBATCH -p hpg-default
#SBATCH -J Nb3Mn8
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --ntasks-per-socket=8
#SBATCH -c 1
#SBATCH --mem=128Gb
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shuan.liu@northeastern.edu

source ~/.bash_aliases; keep_log

module load apptainer/1.3.1

apptainer exec --no-mount /red \\
    --env "LD_LIBRARY_PATH=/blue/m2qm-efrc/shuan.liu.neu/apps/orca/6.1.1/linux_x86-64_shared_openmpi418/lib" \\
    /blue/m2qm-efrc/shuan.liu.neu/containers/orca_3.0.sif \\
    /blue/m2qm-efrc/shuan.liu.neu/apps/orca/6.1.1/linux_x86-64_shared_openmpi418/orca orca.inp 2>&1 > orca.out

rm -rf core.*
"""

job_script = """#!/bin/bash

#SBATCH -A m3346
#SBATCH -J Ce24Cu
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH -t 08:30:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -c 8
#SBATCH --output=output
#SBATCH --error=error
#SBATCH --mail-type=None
#SBATCH --mail-user=shuan.liu@northeastern.edu
#SBATCH --image=shlliu/orca:3.0
#SBATCH --volume=/global/cfs/cdirs/m3346/shlufl/apps/orca/6.1.1/linux_x86-64_shared_openmpi418:/opt/orca

source ~/.bash_aliases; keep_log

shifter --env=LD_LIBRARY_PATH=/opt/orca/lib /opt/orca/orca orca.inp > orca.out

#rm -rf orca.densities orca.gbw *.tmp* # NEVPT2
"""

