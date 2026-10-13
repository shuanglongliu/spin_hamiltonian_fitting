import os
import subprocess
from data import *

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

def get_inp(directories, magmoms):
    ndir = len(directories)
    for i in range(ndir):
        os.chdir(directories[i])
        os.system("pwd")
        magmom = dir2magmom(directories[i], magmoms)
        with open("INCAR", "w") as f:
            f.write(incar.format(magmom=magmom, U=0.0, J=0.0))
        with open("POSCAR", "w") as f:
            f.write(poscar)
        with open("KPOINTS", "w") as f:
            f.write(kpoints)
        with open("vasp.job", "w") as f:
            f.write(job_script)
        os.system("ln -sf " + rootdir + "POTCAR .")
        os.chdir(rootdir)

def submit_job(directories):
    ndir = len(directories)
    for i in range(ndir):
        os.chdir(directories[i])
        os.system("pwd")
        os.system("sbatch vasp.job; sleep 3")
        os.chdir(rootdir)

def get_energies(directories, de0=1.e-8):
    ndir = len(directories)
    nspin = len(directories[0])

    energies = []
    for idir in range(ndir):
        confe = []
        for i, j in enumerate(directories[idir]):
            if j == "d":
                confe.append( magmoms[i]/2)
            else:
                confe.append(-magmoms[i]/2)

        energy = 0
        os.chdir(directories[idir])
        os.system("pwd")
        try:
            out = subprocess.run(["grep", ":", "output"], capture_output=True)
            out = out.stdout.decode("utf-8").split('\n')
            de = out[-2].split()[3]
            print(de)
            if abs(float(de)) < de0:
                try:
                    out = subprocess.run(["grep", "sigma", "OUTCAR"], capture_output=True)
                    energy = float( out.stdout.decode("utf-8").split()[-1] )
                except:
                    pass
        except:
            pass
        os.chdir(rootdir)
        confe.append(energy)

        energies.append(confe)

    with open("energies.dat", "w") as f:
        for idir in range(ndir):
            f.write((nspin*" {:4.1f}" + "{:15.8f}\n").format(*energies[idir]))

    return energies

if __name__ == "__main__":
    directories = ["uuuuuuu", "duuuuuu", "uduuuuu"]
    magmoms = [3, 4, 4, 4, 4, 4, 4]

    #create_dir(directories)
    #get_inp(directories, magmoms)
    #submit_job(directories)
    get_energies(directories, de0=1.e-8)

