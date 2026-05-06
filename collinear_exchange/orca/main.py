"""
Broken-symmetry DFT job manager for the Nd3Mn8 cluster.

Provides the BsDft class for setting up, submitting, backing up, and
resubmitting ORCA broken-symmetry calculations on an HPC cluster via SLURM.

Typical workflow (uncomment the desired block in __main__):
    1. create_work_dir() – create per-config subdirectory
    2. get_input_file()  – write orca.inp / orca.job and set up symlinks
    3. submit_job()      – sbatch orca.job
    4. check_convergence() – inspect orca.out; calls backup() if not converged
    5. resubmit_job()    – resubmit jobs that were backed up but not converged
"""

import os
import subprocess
import csv

from data import (
    root_dir, 
    hs_dir, 
    spin_only_Szs, 
    expected_magmoms,
    effective_spins,
    configs, 
    n_cpus,
    mem_per_core,
    input_file_hs,
    input_file_bs,
    job_script,
)


# Class for one broken-symmetry DFT calculation
class BsDft:
    """Manage a single broken-symmetry DFT calculation for one spin configuration."""

    def __init__(self, config="ud", spin_only_Szs=(0.5, 0.5), restart=False):
        """
        Parameters
        ----------
        config : str
            Spin configuration string ('u'/'d' per magnetic site).
            Must contain at least one 'd' to define a broken-symmetry state.
        spin_only_Szs : sequence of float
            Per-site |ms| values; length must equal len(config).
        input_file : str
            ORCA input file template (may contain {flip_spin} and {final_ms}
            format placeholders).
        job_script : str
            SLURM job script content written verbatim to orca.job.
        """
        global root_dir, input_file, job_script
        # if "d" not in config:
        #     raise ValueError(
        #         "Invalid configuration: 'd' must be present in the configuration string."
        #     )
        self.config = config
        self.spin_only_Szs = spin_only_Szs
        self.restart = restart
        self.work_dir = os.path.join(root_dir, config)
        self.n_sites = len(config)
        self._set_input_file()
        self.job_script = job_script
        self._get_settings()
        self._get_reference_gbw()

    # ─────────────────────────────────────────────────────────────────────────
    # Internal helpers called by __init__
    # ─────────────────────────────────────────────────────────────────────────

    def _get_settings(self):
        """Compute all spin-flip-related parameters for the ORCA input."""
        if self.restart:
            # For the restart input, we want to comment out the FlipSpin and FinalMs directives, since ORCA will read the spin configuration from the broken-symmetry GBW file and applying FlipSpin or FinalMs on top of that will lead to an incorrect state.
            self.flip_spin_tag = "# FlipSpin"
            self.final_ms_tag = "# FinalMs"
        else:
            # For the initial input, we want to specify the FlipSpin directive to define the broken-symmetry state, and FinalMs to reflect that state.
            self.flip_spin_tag = "FlipSpin"
            self.final_ms_tag = "FinalMs"

        # If the number of down spins is bigger than the number of up spins,
        # flip u to d, and d to u, so that the total Ms is always positive. 
        if self.config.count("d") > self.config.count("u"):
            self.config = self.config.replace("u", "x").replace("d", "u").replace("x", "d")

        indices = [i for i, c in enumerate(self.config) if c == "d"]
        self.flip_spin = ", ".join(str(i) for i in indices)

        signed_ms = [
            self.spin_only_Szs[i] if c == "u" else -self.spin_only_Szs[i]
            for i, c in enumerate(self.config)
        ]
        self.final_ms = sum(signed_ms)

        if self.restart:
            S = abs(self.final_ms)
            self.multiplicity = int(2 * S + 1)
        else:
            S = sum(self.spin_only_Szs)
            self.multiplicity = int(2 * S + 1)

    def _get_reference_gbw(self):
        """
        Determine the reference GBW file for the input template based on restart status.
           case 1: initial setup
           case 2: restart after backup
        """
        # Set default reference GBW to the high-spin state
        self.ref_gbw_file = os.path.join(hs_dir, "orca.gbw")
        # For the restart input, we want to use the GBW from the last run of this config
        if self.restart:
            # Find the latest run_N.gbw file in the work_dir
            run_files = [f for f in os.listdir(self.work_dir)
                         if f.startswith("run_") and f.endswith(".gbw")]
            if run_files:
                run_numbers = [int(f.split("_")[1].split(".")[0]) for f in run_files]
                latest_run = max(run_numbers)
                self.ref_gbw_file = os.path.join(self.work_dir, f"run_{latest_run}.gbw")

    def _set_input_file(self):
        """Set input_file to input_file_hs is the state is high-spin,
        otherwise set it to input_file_bs."""
        global input_file_hs, input_file_bs
        if all(c == "u" for c in self.config):
            self.input_file = input_file_hs
        else:
            self.input_file = input_file_bs

    # ─────────────────────────────────────────────────────────────────────────
    # Public workflow methods
    # ─────────────────────────────────────────────────────────────────────────

    def create_work_dir(self):
        """Create the working directory for this configuration if it does not exist."""
        os.makedirs(self.work_dir, exist_ok=True)

    def get_input_file(self):
        """Write orca.inp and orca.job, and create required symlinks in work_dir."""
        wd = self.work_dir; print(wd)

        # Write the ORCA input with the spin-specific placeholders filled in
        with open(os.path.join(wd, "orca.inp"), "w") as f:
            f.write(self.input_file.format(
                n_cpus=n_cpus,
                mem_per_core=mem_per_core,
                flip_spin_tag=self.flip_spin_tag,
                flip_spin=self.flip_spin,
                final_ms_tag=self.final_ms_tag,
                final_ms=self.final_ms,
                multiplicity=self.multiplicity,
            ))

        # Write the SLURM job script verbatim
        with open(os.path.join(wd, "orca.job"), "w") as f:
            f.write(self.job_script)

        # Symlink the shared geometry file (target is relative to the symlink location)
        subprocess.run(["ln", "-sf", "../orca.xyz", os.path.join(wd, "orca.xyz")])

        # Symlink the reference GBW file as the MO starting guess
        subprocess.run(["ln", "-sf", self.ref_gbw_file, os.path.join(wd, "dft.gbw")])

    def submit_job(self):
        """Submit orca.job to SLURM from the working directory."""
        print(self.work_dir)
        subprocess.run(["sbatch", "orca.job"], cwd=self.work_dir)

    def check_convergence(self, verbose=True):
        """
        Check whether the SCF converged.

        Case 1:
            Restart without spin flip: "SCF CONVERGED" should appear once in orca.out.
        Case 2:
            Initial run with spin flip: "SCF CONVERGED" should appear twice in orca.out
            (once for HS calculation, once for BS calculation after the spin flip).
        """
        # if verbose:
        #     print(self.work_dir)

        self.converged = False

        out_path = os.path.join(self.work_dir, "orca.out")
        if os.path.exists(out_path):
            with open(out_path) as f:
                lines = f.readlines()

            # Indicator for initial spin-flip calculation: "F L I P    S P I N"
            if "F L I P    S P I N" in "".join(lines):
                expected_convergence_count = 2
            else:
                expected_convergence_count = 1
    
            convergence_count = sum(1 for line in lines if "SCF CONVERGED" in line)
            if convergence_count == expected_convergence_count:
                self.converged = True

        if verbose:
            if self.converged:
                print(f"Configuration {self.config} is converged.")
            else:
                print(f"Configuration {self.config} is not converged.")

    def backup_files(self):
        """
        Archive completed (or interrupted) run files and prepare for restart.

        Actions:
          - Delete scratch files (.tmp, .bas variants).
          - Copy orca.inp to run_N.inp.
          - Copy orca.out to run_N.out.
          - Move orca.gbw to run_N.gbw.
          - Link run_N.gbw to dft.gbw.
        """
        # Backup only if this run is not converged
        self.check_convergence(verbose=False)
        if self.converged:
            print(f"Configuration {self.config}: converged, skipping backup ...")
            return

        # Remove ORCA scratch files (e.g. *.tmp, *.tmp.1, *.bas)
        wd = self.work_dir
        for fname in os.listdir(wd):
            if ".tmp" in fname or ".bas" in fname:
                os.remove(os.path.join(wd, fname))

        # Proceed only if orca.gbw exists 
        gbw_path = os.path.join(wd, "orca.gbw")
        if not os.path.exists(gbw_path):
            return

        # Check if run_*.gbw files already exist
        run_files = [f for f in os.listdir(wd) if f.startswith("run_") and f.endswith(".gbw")]
        if not run_files:
            # First backup: archive as run_1.*
            run_index = 1
        else:
            # Subsequent backup: increment the run index
            run_numbers = [int(f.split("_")[1].split(".")[0]) for f in run_files]
            run_index = max(run_numbers) + 1
        subprocess.run(["cp", os.path.join(wd, "orca.inp"), os.path.join(wd, f"run_{run_index}.inp")])
        subprocess.run(["cp", os.path.join(wd, "orca.out"), os.path.join(wd, f"run_{run_index}.out")])
        os.rename(gbw_path, os.path.join(wd, f"run_{run_index}.gbw"))
        subprocess.run(["ln", "-sf", f"run_{run_index}.gbw", os.path.join(wd, "dft.gbw")])
        print(f"Configuration {self.config}: backed up the .inp, .out, and .gbw files to run_{run_index}.* and updated dft.gbw symlink.")

    def resubmit_job(self):
        """
        Resubmit a non-converged job.
        """
        # Check convergence
        self.check_convergence(verbose=False)

        # Resubmit only if not converged
        if self.converged:
            print(f"Configuration {self.config}: converged, skipping resubmission ...")
        else:
            print(f"Configuration {self.config}: not converged, resubmitting job ...")
            subprocess.run(["sbatch", "orca.job"], cwd=self.work_dir)

    def check_magnetic_moments(self, moments, magmom_tolerance=0.3):
        """
        Check whether the magnetic moments in the output match the expected configuration.
        Check both the sign (positive for 'u', negative for 'd') and the magnitude (roughly matching the reference expected_magmoms) of the spin populations for each site.
        """
        # Reference magnetic moments 
        global expected_magmoms

        # Are magnetic moments correct?
        self.correct_magmoms = True

        for moment in moments:
            site_index = moment["atomic_index"]
            expected_sign = 1 if self.config[site_index] == "u" else -1
            if moment["spin_population"] * expected_sign < 0:
                self.correct_magmoms = False
                print(f"Configuration {self.config}: site {site_index} (0-based, {moment['element']}) has unexpected spin population in the sign {moment['spin_population']}.")
            expected_mm = expected_magmoms[site_index]
            if abs(abs(moment["spin_population"]) - expected_mm) > magmom_tolerance:
                self.correct_magmoms = False
                print(f"Configuration {self.config}: site {site_index} (0-based, {moment['element']}) has unexpected spin population in the size {moment['spin_population']} (tolerance: {magmom_tolerance} μB).")

    def get_magnetic_moments(self, magmom_tolerance=0.3):
        """
        Extract the final magnetic moments from the orca.out file.
        """

        # Check convergence
        self.check_convergence(verbose=False)
        if not self.converged:
            print(f"Configuration {self.config} is not converged. Stopping ...")
            return

        out_path = os.path.join(self.work_dir, "orca.out")
        with open(out_path) as f:
            lines = f.readlines()

        # Get the magnetic moments from the "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS" section of the output
        moments = []
        sites_parsed = 0
        in_section = False
        for line in lines:
            if "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS" in line:
                moments = [] # Reset moments dict since there are multiple sections 
                sites_parsed = 0
                in_section = True
                continue
            if in_section:
                if line.strip() == "" or line.startswith("----"):
                    continue
                parts = line.split()
                if sites_parsed < self.n_sites:
                    atom_index = int(parts[0])
                    atom_symbol = parts[1].rstrip(":")
                    spin_population = float(parts[3])
                    moments.append({
                        "atomic_index": int(parts[0]),
                        "element": parts[1].rstrip(":"),
                        "spin_population": round(float(parts[3]), 2),
                    })
                    sites_parsed += 1
                else:
                    in_section = False
                    # Do not break here since there may be multiple sections in the output (e.g.
                    # for the initial HS calculation and the subsequent BS calculation), and we
                    # want to capture the moments from the final BS section.

        # Check the extracted magnetic moments against the expected configuration
        self.check_magnetic_moments(moments, magmom_tolerance=magmom_tolerance)

        # Save the moments to a .csv file in the work_dir
        with open(os.path.join(self.work_dir, "magnetic_moments.csv"), "w", newline="") as csvfile:
            fieldnames = ["atomic_index", "element", "spin_population"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for moment in moments:
                writer.writerow(moment)

    def get_energy(self, magmom_tolerance=0.3):
        """
        Extract the final SCF energy from the orca.out file.
        """
        self.energy = 0.0

        # Check convergence
        self.check_convergence(verbose=False)
        if not self.converged:
            return None

        # Check magnetic moments
        self.get_magnetic_moments(magmom_tolerance=magmom_tolerance)
        if not self.correct_magmoms:
            return None

        out_path = os.path.join(self.work_dir, "orca.out")
        with open(out_path) as f:
            lines = f.readlines()

        for line in reversed(lines):
            if "FINAL SINGLE POINT ENERGY" in line:
                self.energy = float(line.split()[-1])
                break

# Class for multiple broken-symmetry DFT calculations
class BsDftManager:
    """Manage multiple broken-symmetry DFT calculations across different configurations."""

    def __init__(self, configs, spin_only_Szs, restart=False):
        self.configs = configs
        self.spin_only_Szs = spin_only_Szs
        self.restart = restart
        self.jobs = [BsDft(config=config, spin_only_Szs=spin_only_Szs, restart=restart) for config in configs]

    def create_work_dirs(self):
        for job in self.jobs:
            job.create_work_dir()

    def get_input_files(self):
        for job in self.jobs:
            job.get_input_file()

    def submit_jobs(self):
        for job in self.jobs:
            job.submit_job()

    def check_convergence(self):
        for job in self.jobs:
            job.check_convergence()

    def backup_files(self):
        for job in self.jobs:
            job.backup_files()

    def resubmit_jobs(self):
        for job in self.jobs:
            job.resubmit_job()

    def get_magnetic_moments(self, magmom_tolerance=0.3):
        for job in self.jobs:
            job.get_magnetic_moments(magmom_tolerance=magmom_tolerance)

    def get_energies(self, magmom_tolerance=0.3):
        global effective_spins
        # Save energies to a data file energies.dat 
        with open(os.path.join(root_dir, "energies.dat"), "w") as f:
            # Add a header line with the site labels and the energy label
            site_labels = "   ".join(f"S{i+1}" for i in range(len(self.spin_only_Szs)))
            f.write(f"# {site_labels}   Energy (Ha)\n")
            for job in self.jobs:
                job.get_energy(magmom_tolerance=magmom_tolerance)
                if job.energy is not None:
                    config_str = " ".join(f"{+effective_spins[i]:5.1f}" if c == "u" else f"{-effective_spins[i]:5.1f}" for i, c in enumerate(job.config))
                    f.write(f"{config_str} {job.energy:.12f}\n")
        print(f"Saved energies for all configurations to {os.path.join(root_dir, 'energies.dat')}")

if __name__ == "__main__":

    # Instantiate the manager for all configurations
    manager = BsDftManager(configs=configs, spin_only_Szs=spin_only_Szs, restart=False)

    # Step 1: Create work directories for all configurations
    # manager.create_work_dirs()

    # Step 2: Write input files and set up symlinks
    manager.get_input_files()

    # Step 3: Submit jobs to SLURM
    # manager.submit_jobs()

    # Step 4: Check convergence and back up files for non-converged jobs
    # manager.check_convergence()
    # manager.backup_files()

    # Step 5: Resubmit non-converged jobs after backup
    # manager.resubmit_jobs()

    # Step 6: After convergence, extract magnetic moments and energies
    # manager.get_magnetic_moments(magmom_tolerance=0.5)
    # manager.get_energies(magmom_tolerance=0.3)

