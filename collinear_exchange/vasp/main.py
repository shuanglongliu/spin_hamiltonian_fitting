import os
import subprocess
from data import root_dir, incar, poscar, kpoints, job_script

def create_dir(directories):
    ndir = len(directories)
    for i in range(ndir):
        os.system("mkdir " + directories[i])
    return

def dir2magmom(dirname, magmoms):
    magmoms_str = ""
    for i, j in enumerate(dirname):
        if j == "d":
            magmoms_str = magmoms_str + " {:4.1f}".format(-magmoms[i])
        else:
            magmoms_str = magmoms_str + " {:4.1f}".format( magmoms[i])
    magmoms_str = magmoms_str + " 1000*0.0"
    return magmoms_str

def get_total_magmom(dirname, magmoms):
    total_magmom = 0
    for i, j in enumerate(dirname):
        if j == "d":
            total_magmom = total_magmom - magmoms[i]
        else:
            total_magmom = total_magmom + magmoms[i]
    return total_magmom

def get_inp(directories, magmoms):
    ndir = len(directories)
    for i in range(ndir):
        os.chdir(directories[i])
        os.system("pwd")
        magmom = dir2magmom(directories[i], magmoms)
        nupdown = get_total_magmom(directories[i], magmoms)
        with open("INCAR", "w") as f:
            f.write(incar.format(magmom=magmom, nupdown=nupdown, U=0.0, J=0.0))
        with open("POSCAR", "w") as f:
            f.write(poscar)
        with open("KPOINTS", "w") as f:
            f.write(kpoints)
        with open("vasp.job", "w") as f:
            f.write(job_script)
        os.system("ln -sf " + os.path.join(root_dir, "POTCAR") + " .")
        os.chdir(root_dir)

def submit_job(directories):
    ndir = len(directories)
    for i in range(ndir):
        os.chdir(directories[i])
        os.system("pwd")
        os.system("sbatch vasp.job")
        os.chdir(root_dir)

def check_convergences(directories):
    ndir = len(directories)

    for idir in range(ndir):
        os.chdir(directories[idir])
        os.system("pwd")
        try:
            out = subprocess.run(["grep", "EDIFF is reached", "OUTCAR"], text=True, capture_output=True)
            print(out.stdout)
            if out.returncode == 0:
                print("Convergence reached in " + directories[idir])
            else:
                print("Convergence NOT reached in " + directories[idir])
        except:
            pass
        os.chdir(root_dir)

    return 

def check_magmoms(directories):
    ndir = len(directories)
    n_sites = len(directories[0])

    for idir in range(ndir):
        os.chdir(directories[idir])
        os.system("pwd")
        try:
            with open("OUTCAR", "r") as f:
                lines = f.readlines()
            magmoms = []
            for i, line in enumerate(lines):
                if "magnetization (x)" in line:
                    magmoms = [float(lines[i+j].split()[-1]) for j in range(4, 4+n_sites)]
            print(magmoms)
        except:
            pass
        os.chdir(root_dir)

    return 

def get_energies(directories):
    global effective_spins

    ndir = len(directories)
    nspin = len(directories[0])

    energies = []
    for idir in range(ndir):
        confe = []
        for i, j in enumerate(directories[idir]):
            if j == "d":
                confe.append(-effective_spins[i])
            else:
                confe.append( effective_spins[i])

        energy = 0
        os.chdir(directories[idir])
        os.system("pwd")
        try:
            out = subprocess.run(["grep", "sigma", "OUTCAR"], capture_output=True)
            energy = float( out.stdout.decode("utf-8").split()[-1] )
        except:
            pass
        os.chdir(root_dir)
        confe.append(energy)

        energies.append(confe)

    with open("energies.dat", "w") as f:
        # Write the header
        f.write("# ")
        for i in range(nspin):
            f.write("S{:d}   ".format(i+1))
        f.write("Energy (eV)\n")
        for idir in range(ndir):
            f.write((nspin*" {:5.1f}" + " {:15.8f}\n").format(*energies[idir]))

    return energies

if __name__ == "__main__":
    # There are 12 exchange couplings. At least 13 magnetic configurations are needed to determine them.
    directories = ["uuuuuuu", "uduuuuu", "uuduuuu", "uuuduuu", "udduuuu", "ududuuu", "uudduuu", "uuuuduu", "uuuuudu", "uuuuuud", "uuuuddu", "uuuudud", "uuuuudd"]
    magmoms = [2, 2, 2, 2, 2, 2, 2]
    effective_spins = [1, 1, 1, 1, 1, 1, 1]

    # create_dir(directories)
    # get_inp(directories, magmoms)
    # submit_job(directories)
    # check_convergences(directories)
    # check_magmoms(directories)
    get_energies(directories)

