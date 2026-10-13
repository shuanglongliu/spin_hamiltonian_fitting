"""
Fit the coefficients of the cross- and dot-product spin terms in the electric polarization.

    p = sum_{ij} [ alpha_ij (Si (Sj . eij) - (Si . eij) Sj) + pi_ij (Si . Sj) ]

eij is the unit vector from site i to site j. The fit minimizes |A x - y|^2,
subject to C x = z when constraints are given, where y holds the changes in
the dipole of every configuration relative to the first one and

    x = [alpha_pair1, ..., alpha_pairn,
         pi_pair1_x, pi_pair1_y, pi_pair1_z, ..., pi_pairn_x, pi_pairn_y, pi_pairn_z]

Inputs (in the working directory): input.yaml and dipoles.dat
(columns: theta, phi per site in degrees, then px, py, pz in e Ang).
Outputs: log and fitting_quality.pdf.

Usage: python main.py [working directory]
"""

import os
import sys

import numpy as np
import yaml
from scipy.linalg import null_space
from scipy.spatial.transform import Rotation


def sph2cart(r: float, theta: float, phi: float) -> np.ndarray:
    """Spherical (angles in degrees) to Cartesian coordinates."""
    theta, phi = np.deg2rad(theta), np.deg2rad(phi)
    return r * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def euler_zyz(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Rotation matrix Rz(gamma) Ry(beta) Rz(alpha) (extrinsic zyz, degrees)."""
    return Rotation.from_euler("zyz", [alpha, beta, gamma], degrees=True).as_matrix()


def read_yaml() -> dict:
    """Read spins, pairs, bond unit vectors, and constraints from input.yaml."""
    with open("input.yaml") as f:
        data = yaml.safe_load(f)
    positions = np.array(data["positions"], dtype=float)
    pairs = data["pairs"]
    eij = np.array([positions[j - 1] - positions[i - 1] for i, j in pairs])
    eij /= np.linalg.norm(eij, axis=1)[:, None]
    return {
        "spins": np.array(data["spins"], dtype=float),
        "pairs": pairs,
        "eij": eij,
        "constraints_alpha": data.get("constraints_alpha") or [],
        "constraints_pi": data.get("constraints_pi") or [],
    }


def read_dipoles(spins: np.ndarray):
    """Read dipoles.dat.

    Returns
    -------
    S : np.ndarray
        Spin vectors, shape (n_conf, n_site, 3).
    P : np.ndarray
        Dipoles, shape (n_conf, 3).
    """
    data = np.loadtxt("dipoles.dat", ndmin=2)
    n_site = len(spins)
    S = np.array([[sph2cart(spins[j], row[2 * j], row[2 * j + 1]) for j in range(n_site)] for row in data])
    P = data[:, 2 * n_site:2 * n_site + 3]
    return S, P


def get_y(P: np.ndarray, i_conf_ref: int = 0) -> np.ndarray:
    """Changes in the dipole relative to the reference configuration, flattened (conf, xyz)."""
    return np.delete(P - P[i_conf_ref], i_conf_ref, axis=0).flatten()


def get_A(S: np.ndarray, eij: np.ndarray, pairs: list, i_conf_ref: int = 0) -> np.ndarray:
    """Coefficient matrix A, shape (3 * (n_conf - 1), 4 * n_pair)."""
    n_conf = S.shape[0]
    n_pair = len(pairs)

    # B holds the coefficients of every configuration before subtracting the reference
    B = np.zeros((n_conf, 3, 4 * n_pair))
    for k, (i, j) in enumerate(pairs):
        Si, Sj = S[:, i - 1], S[:, j - 1]
        B[:, :, k] = Si * (Sj @ eij[k])[:, None] - (Si @ eij[k])[:, None] * Sj
        SiSj = np.einsum("ca,ca->c", Si, Sj)
        for a in range(3):
            B[:, a, n_pair + 3 * k + a] = SiSj

    return np.delete(B - B[i_conf_ref], i_conf_ref, axis=0).reshape(-1, 4 * n_pair)


def get_constraints(constraints_alpha: list, constraints_pi: list, n_pair: int):
    """Constraint matrix C and right-hand side z.

    constraints_alpha: [pair_i, pair_j] means alpha(pair_i) = alpha(pair_j).
    constraints_pi: [pair_i, pair_j, alpha, beta, gamma] means pi(pair_j) = R pi(pair_i),
    with R from extrinsic zyz Euler angles in degrees.
    """
    n_c = len(constraints_alpha) + 3 * len(constraints_pi)
    C = np.zeros((n_c, 4 * n_pair))

    for row, (i_pair, j_pair) in enumerate(constraints_alpha):
        C[row, i_pair - 1] = 1.0
        C[row, j_pair - 1] = -1.0

    row = len(constraints_alpha)
    for i_pair, j_pair, alpha, beta, gamma in constraints_pi:
        rot = euler_zyz(alpha, beta, gamma)
        for a in range(3):
            C[row, n_pair + 3 * (i_pair - 1):n_pair + 3 * i_pair] = rot[a]
            C[row, n_pair + 3 * (j_pair - 1) + a] = -1.0
            row += 1

    return C, np.zeros(n_c)


def solve(A: np.ndarray, y: np.ndarray, C: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Minimize |A x - y| subject to C x = z.

    x is written as x_p + N w, where x_p solves C x = z and the columns of N span
    the null space of C, and w is found by least squares. If the data do not
    determine every free parameter, the minimum-norm solution is returned and a
    warning is printed.
    """
    if C.shape[0] == 0:
        x_p = np.zeros(A.shape[1])
        N = np.eye(A.shape[1])
    else:
        x_p = np.linalg.lstsq(C, z, rcond=None)[0]
        N = null_space(C)
    w, _, rank, _ = np.linalg.lstsq(A @ N, y - A @ x_p, rcond=None)
    if rank < N.shape[1]:
        print(f"Warning: {N.shape[1] - rank} combination(s) of the {N.shape[1]} free parameters "
              f"are not determined by the data; returning the minimum-norm solution.")
    return x_p + N @ w


def print_results(x: np.ndarray, pairs: list):
    n_pair = len(pairs)
    with open("log", "w") as f:
        for k, (i, j) in enumerate(pairs):
            f.write(f"Pair {i:3d} {j:3d}\n\n")
            f.write(f"alpha = {x[k]:12.6f} e Ang\n")
            pi = x[n_pair + 3 * k:n_pair + 3 * k + 3]
            f.write(f"pi = ({pi[0]:12.6f} {pi[1]:12.6f} {pi[2]:12.6f}) e Ang\n")
            f.write("\n")


def analysis(x: np.ndarray, y: np.ndarray, A: np.ndarray):
    """Append the deviations of the fit to log and plot DFT vs fit."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y_fit = A @ x
    dy = np.abs(y_fit - y)
    with open("log", "a") as f:
        f.write(f"\n\nMaximal deviation: {dy.max():12.6f} e Ang.\n")
        f.write(f"Average deviation: {dy.mean():12.6f} e Ang.\n")

    index = np.arange(1, len(y) + 1)
    fig, ax = plt.subplots()
    ax.scatter(index, y, label="DFT")
    ax.scatter(index, y_fit, label="Fit")
    ax.set_xlabel(r"Index of $\delta P_\alpha$ for all spin configurations")
    ax.set_ylabel(r"$\delta P_\alpha$ (e $\AA$)")
    ax.legend(loc="lower right", frameon=False)
    fig.savefig("fitting_quality.pdf")
    plt.close(fig)


def main():
    inp = read_yaml()
    S, P = read_dipoles(inp["spins"])

    y = get_y(P)
    A = get_A(S, inp["eij"], inp["pairs"])
    C, z = get_constraints(inp["constraints_alpha"], inp["constraints_pi"], len(inp["pairs"]))
    x = solve(A, y, C, z)

    print_results(x, inp["pairs"])
    analysis(x, y, A)
    return x


if __name__ == "__main__":
    if len(sys.argv) > 1:
        os.chdir(sys.argv[1])
    main()
