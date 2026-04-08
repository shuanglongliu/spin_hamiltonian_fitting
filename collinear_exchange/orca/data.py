"""
Static data for broken-symmetry DFT calculations on the Nd3Mn8 cluster.

Defines spin configurations, per-site spin quantum numbers, ORCA input
file templates, and the SLURM job script used by main.py.
"""

import os

# Paths anchored to this file's location
root_dir = os.path.dirname(os.path.abspath(__file__))

# High-spin (all-up) reference calculation used as MO starting guess
hs_dir = os.path.join(root_dir, "uuuuuuuuuuu")

# Per-site spin quantum numbers |ms| for the 11 magnetic centers:
spin_only_Szs = [1.5, 1.0, 1.5, 2, 2, 2, 2, 2, 2, 2, 2]

# Effective spin quantum numbers
effective_spins = [4.5, 4.0, 4.5, 2, 2, 2, 2, 2, 2, 2, 2]

# Expected local magnetic moments (in μB) for the HS reference, from Loewdin analysis:
expected_magmoms = [
3.046566,
2.827955,
3.048257,
3.777876,
3.662067,
3.776683,
3.719454,
3.712534,
3.782360,
3.662376,
3.765104 
]

# Broken-symmetry spin configurations for the 11-site cluster.
# Each character gives the local spin orientation at that site:
#   'u' = spin-up (ms = +|ms|),  'd' = spin-down (ms = −|ms|)
# The first block contains coupled/mixed states; the remaining blocks
# enumerate all single-site and two-site spin flips from the HS reference.
configs = [
    "uuuudduuddu",
    "uuuuudduudd",
    "uuudddduuuu",
    "uuudduuuudd",
    "uuuduududdu",
    "uuududududu",
    "uuududuudud",
    "ddduuuuuudu",
    "uuudduuuuud",
    "ddduuduuuud",
    # Single-site spin-down flips (one 'd' at each position)
    "duuuuuuuuuu",
    "uduuuuuuuuu",
    "uuduuuuuuuu",
    "uuuduuuuuuu",
    "uuuuduuuuuu",
    "uuuuuduuuuu",
    "uuuuuuduuuu",
    "uuuuuuuduuu",
    "uuuuuuuuduu",
    "uuuuuuuuudu",
    "uuuuuuuuuud",
    # Two-site spin-down flips (first site fixed as 'd', second site varies)
    "dduuuuuuuuu",
    "duduuuuuuuu",
    "duuduuuuuuu",
    "duuuduuuuuu",
    "duuuuduuuuu",
    "duuuuuduuuu",
    "duuuuuuduuu",
    "duuuuuuuduu",
    "duuuuuuuudu",
    "duuuuuuuuud",
    "udduuuuuuuu",
    "ududuuuuuuu",
    "uduuduuuuuu",
    "uduuuduuuuu",
    "uduuuuduuuu",
    "uduuuuuduuu",
    "uduuuuuuduu",
    "uduuuuuuudu",
    "uduuuuuuuud",
    "uudduuuuuuu",
    "uududuuuuuu",
    "uuduuduuuuu",
    "uuduuuduuuu",
    "uuduuuuduuu",
    "uuduuuuuduu",
    "uuduuuuuudu",
    "uuduuuuuuud",
    "uuudduuuuuu",
    "uuududuuuuu",
    "uuuduuduuuu",
    "uuuduuuduuu",
    "uuuduuuuduu",
    "uuuduuuuudu",
    "uuuduuuuuud",
    "uuuudduuuuu",
    "uuuududuuuu",
    "uuuuduuduuu",
    "uuuuduuuduu",
    "uuuuduuuudu",
    "uuuuduuuuud",
    "uuuuudduuuu",
    "uuuuududuuu",
    "uuuuuduuduu",
    "uuuuuduuudu",
    "uuuuuduuuud",
    "uuuuuudduuu",
    "uuuuuududuu",
    "uuuuuuduudu",
    "uuuuuuduuud",
    "uuuuuuudduu",
    "uuuuuuududu",
    "uuuuuuuduud",
    "uuuuuuuuddu",
    "uuuuuuuudud",
    "uuuuuuuuudd",
]

# Push the HS reference configuration to the front of the list for convenience
configs.insert(0, "uuuuuuuuuuu")

n_configs = len(configs)

# ─────────────────────────────────────────────────────────────────────────────
# ORCA input file templates
# ─────────────────────────────────────────────────────────────────────────────

input_file = """
# --- Method & relativistic treatment ---
# ZORA scalar relativistics are essential for lanthanides
# TPSSh gives reliable J values for transition/lanthanide systems
! TPSSh ZORA DefGrid3 TightSCF KDIIS MORead # SlowConv NoTRAH 

%basis
    # All-electron ZORA basis for lanthanide atoms — must use SARC
    NewGTO Nd "SARC-ZORA-TZVP" end

    # ZORA-recontracted triple-zeta for transition metal and light atoms
    NewGTO Mn "ZORA-def2-TZVP"  end
    NewGTO C  "ZORA-def2-TZVP"  end
    NewGTO N  "ZORA-def2-TZVP"  end
    NewGTO O  "ZORA-def2-TZVP"  end
    NewGTO H  "ZORA-def2-SVP"   end   # SVP for H is sufficient
end

# --- Parallelization ---
%pal nprocs 16 end

# --- Memory per core in MB ---
%maxcore 8000

# --- SCF options ---
%scf
    MaxIter   500

    # Use larger DIIS to help convergence with open-shell lanthanide
    DIISMaxEq 15

    # Trusted region augmentation Hessian (TRAH): 2nd order SCF
    AutoTRAH  true

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
#  CHARGE and MULTIPLICITY  (adjust to your system!)
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
