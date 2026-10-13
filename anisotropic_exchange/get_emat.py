
import numpy as np
from ase.io import read, write

def print_emat_array(emat):
    print("emats_file[0] = [ \\")
    print("        [ {:15.8f}, {:15.8f}, {:15.8f} ], \\".format(*emat[0]))
    print("        [ {:15.8f}, {:15.8f}, {:15.8f} ], \\".format(*emat[1]))
    print("        [ {:15.8f}, {:15.8f}, {:15.8f} ], \\".format(*emat[2]))
    print("        ]")

def print_emat_vector(emat):
    emat_flat = emat.flatten()
    print("reference_frame: [{:15.8f}, {:15.8f}, {:15.8f}, {:15.8f}, {:15.8f}, {:15.8f}, {:15.8f}, {:15.8f}, {:15.8f}]".format(*emat_flat))

def get_emat_exchange(i = 1, j = 2, form="vector"):

    atoms = read("left_hand_relaxed.vasp")
    positions = atoms.get_positions()
    ex = positions[j-1] - positions[i-1]
    ex /= np.linalg.norm(ex)
    ez = np.array([0., 0., 1.])
    ey = np.cross(ez, ex)
    emat = np.array([ex, ey, ez])

    print("Pair {:3d} {:3d}".format(i, j))
    if form == "vector":
        print_emat_vector(emat)
    else:
        print_emat_array(emat)

if __name__ == "__main__":

    get_emat_exchange(i = 1, j = 2)

