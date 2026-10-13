"""
Fit full 3x3 exchange coupling matrices from noncollinear DFT energies.

Spin Hamiltonian: H = factor * sum_{ij} Si^T J_ij Sj

Each J_ij is expressed in a local reference frame attached to the pair (i, j).
Every energy pair (conf1, conf2) in input.yaml gives one equation

    E2 - E1 = factor * sum_{ij} sum_{ab} J_ij^{ab} (S2_ia S2_jb - S1_ia S1_jb)

with the spin components rotated into the frame of pair (i, j).
The unknowns are ordered pair-1_xx, pair-1_xy, ..., pair-n_zz.

Usage: python main.py [directory containing input.yaml]
"""

import os
import sys

import numpy as np
import yaml


def sph2cart(r: float, theta: float, phi: float) -> np.ndarray:
    """Spherical (angles in degrees) to Cartesian coordinates."""
    theta, phi = np.deg2rad(theta), np.deg2rad(phi)
    return r * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def read_input(dirname: str) -> dict:
    """Read input.yaml.

    Returns
    -------
    dict
        spins : (n_site,) spin quantum numbers
        spin_pairs : list of 1-indexed [i, j] pairs, one per J matrix
        reference_frames : list of (3, 3) arrays whose columns are the local axes ex, ey, ez
        energy_pairs : the "energies" entries of input.yaml
    """
    with open(os.path.join(dirname, "input.yaml")) as f:
        data = yaml.safe_load(f)
    return {
        "spins": np.array(data["spins"], dtype=float),
        "spin_pairs": [list(entry["pair"]) for entry in data["exchange"]],
        # reference_frame lists ex, ey, ez one after another, so the rows of the
        # reshaped array are the local axes; transpose to make them the columns.
        "reference_frames": [np.array(entry["reference_frame"], dtype=float).reshape(3, 3).T
                             for entry in data["exchange"]],
        "energy_pairs": data["energies"],
    }


def get_spin_vectors(conf: list, spins: np.ndarray) -> np.ndarray:
    """Spin vectors (n_site, 3) of one configuration [theta_1, phi_1, ..., energy]."""
    return np.array([sph2cart(S, conf[2 * i], conf[2 * i + 1]) for i, S in enumerate(spins)])


def get_coefficient_matrix(inp: dict, factor: float):
    """Build the coefficient matrix C and the energy differences dE (meV).

    Returns
    -------
    C : np.ndarray
        Shape (n_energy_pair, 9 * n_spin_pair).
    dE : np.ndarray
        E2 - E1 for each energy pair in meV (energies are read in eV).
    """
    spins = inp["spins"]
    spin_pairs = inp["spin_pairs"]
    frames = inp["reference_frames"]
    n_site = len(spins)

    C = np.zeros((len(inp["energy_pairs"]), 9 * len(spin_pairs)))
    dE = np.zeros(len(inp["energy_pairs"]))
    for row, energy_pair in enumerate(inp["energy_pairs"]):
        # The angles are given in the frame of the listed pair
        frame_from = frames[spin_pairs.index(list(energy_pair["pair"]))]
        S1 = get_spin_vectors(energy_pair["conf1"], spins)
        S2 = get_spin_vectors(energy_pair["conf2"], spins)
        for k, (i, j) in enumerate(spin_pairs):
            # Rotate the spins into the frame of pair k
            T = frames[k].T @ frame_from
            S1_k = S1 @ T.T
            S2_k = S2 @ T.T
            C[row, 9 * k:9 * (k + 1)] = factor * (np.outer(S2_k[i - 1], S2_k[j - 1])
                                                  - np.outer(S1_k[i - 1], S1_k[j - 1])).flatten()
        dE[row] = 1000 * (energy_pair["conf2"][2 * n_site] - energy_pair["conf1"][2 * n_site])
    return C, dE


def get_Js(inp: dict, factor: float = -2.0, quantum_correction: bool = False) -> np.ndarray:
    """Solve for the J matrix components (meV).

    An exact solve is used when there are exactly 9 energy pairs per spin pair,
    least squares when there are more.

    Raises
    ------
    ValueError
        If the system is underdetermined.
    """
    C, dE = get_coefficient_matrix(inp, factor)
    n_eq, n_J = C.shape
    if n_eq < n_J:
        raise ValueError(f"Underdetermined system: {n_eq} energy pair(s) for {n_J} unknown(s).")

    if n_eq == n_J:
        Js = np.linalg.solve(C, dE)
    else:
        Js = np.linalg.lstsq(C, dE, rcond=None)[0]
        print(f"Least-squares fit of {n_eq} energy pairs for {n_J} unknowns: "
              f"residual = {np.linalg.norm(C @ Js - dE):.3e} meV")

    if quantum_correction:
        # See Ruiz 2003. This is purely empirical.
        spins = inp["spins"]
        for k, (i, j) in enumerate(inp["spin_pairs"]):
            Si, Sj = spins[i - 1], spins[j - 1]
            Js[9 * k:9 * (k + 1)] *= 2 * Si * Sj / (2 * Si * Sj + min(Si, Sj))
    return Js


def print_Js(Js: np.ndarray, iformat: int = 2):
    """Print the J components.

    iformat = 1: labelled, one block per pair; iformat = 2: all components on one line.
    """
    labels = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]
    if iformat == 1:
        for k in range(len(Js) // 9):
            print(f"Pair {k + 1:3d}")
            for label, J in zip(labels, Js[9 * k:9 * (k + 1)]):
                print(f"J_{label} = {J:9.3f}")
    else:
        print("".join(f" {J:9.3f}" for J in Js))


if __name__ == "__main__":
    dirname = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    inp = read_input(dirname)
    Js = get_Js(inp, factor=-2.0, quantum_correction=False)
    print_Js(Js, iformat=2)
