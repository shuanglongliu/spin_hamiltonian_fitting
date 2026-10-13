"""
Fit the full spin-driven polarization tensors P^{ij} for every pair, with symmetry constraints.

    p_gamma = sum_{ij} sum_{alpha beta} S_{i alpha} P^{ij}_{gamma alpha beta} S_{j beta}

The fit minimizes |A x - y|^2 subject to C x = z, where each constraint
[pair_i, pair_j, alpha, beta, gamma] in input.yaml imposes

    P^{pair_j}_{gamma alpha beta} = sum R_{gamma gamma'} R_{alpha alpha'} R_{beta beta'} P^{pair_i}_{gamma' alpha' beta'}

The elements are ordered with beta fastest, then alpha, then gamma, then the pair:
x[27 * k + 9 * gamma + 3 * alpha + beta].

Inputs (in the working directory): input.yaml and dipoles.dat
(columns: theta, phi per site in degrees, then px, py, pz in e Ang).
Outputs: P_total.log, fit.log, and fitting_quality.pdf.

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
    """Read spins, pairs, and constraints from input.yaml."""
    with open("input.yaml") as f:
        data = yaml.safe_load(f)
    return {
        "spins": np.array(data["spins"], dtype=float),
        "pairs": data["pairs"],
        "constraints": data.get("constraints") or [],
    }


def read_dipoles(spins: np.ndarray):
    """Read dipoles.dat; return spin vectors (n_conf, n_site, 3) and dipoles (n_conf, 3)."""
    data = np.loadtxt("dipoles.dat", ndmin=2)
    n_site = len(spins)
    S = np.array([[sph2cart(spins[j], row[2 * j], row[2 * j + 1]) for j in range(n_site)] for row in data])
    P = data[:, 2 * n_site:2 * n_site + 3]
    return S, P


def get_y(P: np.ndarray, i_conf_ref: int = 0) -> np.ndarray:
    """Changes in the dipole relative to the reference configuration, flattened (conf, gamma)."""
    return np.delete(P - P[i_conf_ref], i_conf_ref, axis=0).flatten()


def get_A(S: np.ndarray, pairs: list, i_conf_ref: int = 0) -> np.ndarray:
    """Coefficient matrix A, shape (3 * (n_conf - 1), 27 * n_pair)."""
    n_conf = S.shape[0]
    B = np.zeros((n_conf, 3, 27 * len(pairs)))
    for k, (i, j) in enumerate(pairs):
        SiSj = np.einsum("ca,cb->cab", S[:, i - 1], S[:, j - 1]).reshape(n_conf, 9)
        for gamma in range(3):
            B[:, gamma, 27 * k + 9 * gamma:27 * k + 9 * (gamma + 1)] = SiSj
    return np.delete(B - B[i_conf_ref], i_conf_ref, axis=0).reshape(-1, 27 * len(pairs))


def get_constraints(constraints: list, n_pair: int):
    """Constraint matrix C (27 rows per constraint) and right-hand side z."""
    C = np.zeros((27 * len(constraints), 27 * n_pair))
    for ic, (i_pair, j_pair, alpha, beta, gamma) in enumerate(constraints):
        rot = euler_zyz(alpha, beta, gamma)
        rows = slice(27 * ic, 27 * (ic + 1))
        C[rows, 27 * (i_pair - 1):27 * i_pair] = np.kron(rot, np.kron(rot, rot))
        C[rows, 27 * (j_pair - 1):27 * j_pair] = -np.eye(27)
    return C, np.zeros(C.shape[0])


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
    """Write the 3x3 matrix P_gamma for every gamma and pair to P_total.log."""
    with open("P_total.log", "w") as f:
        for k, (i, j) in enumerate(pairs):
            f.write(f"Pair {i:3d} {j:3d}\n\n")
            for gamma, direction in enumerate("xyz"):
                f.write(f"The {direction} component of coupling matrix\n\n")
                for alpha in range(3):
                    index = 27 * k + 9 * gamma + 3 * alpha
                    f.write("{:12.6f} {:12.6f} {:12.6f}\n".format(*x[index:index + 3]))
                f.write("\n")
            f.write("\n")


def analysis(x: np.ndarray, y: np.ndarray, A: np.ndarray):
    """Write deviation statistics to fit.log and plot DFT vs fit."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y_fit = A @ x
    dy = np.abs(y_fit - y)
    with open("fit.log", "w") as f:
        f.write(f"Maximal deviation: {dy.max():12.6f} e Ang.\n")
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
    A = get_A(S, inp["pairs"])
    C, z = get_constraints(inp["constraints"], len(inp["pairs"]))
    x = solve(A, y, C, z)

    print_results(x, inp["pairs"])
    analysis(x, y, A)
    return x


if __name__ == "__main__":
    if len(sys.argv) > 1:
        os.chdir(sys.argv[1])
    main()
