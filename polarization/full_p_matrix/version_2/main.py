"""
Fit the full spin-driven polarization tensors P^{ij}, using only the irreducible pairs as unknowns.

    p_gamma = sum_{ij} sum_{alpha beta} S_{i alpha} P^{ij}_{gamma alpha beta} S_{j beta}

For an irreducible pair ij the unknowns are the 27 elements of P^{ij}. A pair ab
related to an irreducible pair ij by a rotation R uses

    P^{ab}_{gamma alpha beta} = sum R_{gamma gamma'} R_{alpha alpha'} R_{beta beta'} P^{ij}_{gamma' alpha' beta'}

The elements are ordered with beta fastest, then alpha, then gamma, then the pair:
x[27 * k + 9 * gamma + 3 * alpha + beta].

Inputs (in the working directory): input.yaml and dipoles.dat
(columns: theta, phi per site in degrees, then px, py, pz in e Ang).

Usage: python main.py [working directory]
"""

import os
import sys

import numpy as np
import yaml
from scipy.spatial.transform import Rotation


def sph2cart(r: float, theta: float, phi: float) -> np.ndarray:
    """Spherical (angles in degrees) to Cartesian coordinates."""
    theta, phi = np.deg2rad(theta), np.deg2rad(phi)
    return r * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def euler_zyz(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Rotation matrix Rz(gamma) Ry(beta) Rz(alpha) (extrinsic zyz, degrees)."""
    return Rotation.from_euler("zyz", [alpha, beta, gamma], degrees=True).as_matrix()


def read_yaml() -> dict:
    """Read input.yaml.

    Returns
    -------
    dict
        spins, pairs, irr_pairs (1-indexed indices into pairs),
        rmats (one (3, 3) array per irreducible pair, columns = local axes ex, ey, ez),
        transformations ([pair_i, pair_j, alpha, beta, gamma]: P(pair_j) = R P(pair_i)).
    """
    with open("input.yaml") as f:
        data = yaml.safe_load(f)
    return {
        "spins": np.array(data["spins"], dtype=float),
        "pairs": data["pairs"],
        "irr_pairs": data["irr_pairs"],
        # emats list ex, ey, ez one after another; rmat is the transpose of emat
        "rmats": [np.array(e, dtype=float).reshape(3, 3).T for e in data["emats"]],
        "transformations": data.get("transformations") or [],
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


def get_A(S: np.ndarray, inp: dict, i_conf_ref: int = 0) -> np.ndarray:
    """Coefficient matrix A, shape (3 * (n_conf - 1), 27 * n_irrp)."""
    n_conf = S.shape[0]
    irr_pairs = inp["irr_pairs"]
    transformations = {t[1]: t for t in inp["transformations"]}

    # B holds the coefficients of every configuration before subtracting the reference
    B = np.zeros((n_conf, 3, 27 * len(irr_pairs)))
    for k_pair, (i_site, j_site) in enumerate(inp["pairs"], start=1):
        Si, Sj = S[:, i_site - 1], S[:, j_site - 1]
        if k_pair in irr_pairs:
            rot = np.eye(3)
            i_irrp = irr_pairs.index(k_pair)
        else:
            i_pair, _, alpha, beta, gamma = transformations[k_pair]
            rot = euler_zyz(alpha, beta, gamma)
            i_irrp = irr_pairs.index(i_pair)
        # B[c, gamma, (gamma', alpha', beta')] = R_{gamma gamma'} (R^T Si)_{alpha'} (R^T Sj)_{beta'}
        block = np.einsum("gh,ca,cb->cghab", rot, Si @ rot, Sj @ rot).reshape(n_conf, 3, 27)
        B[:, :, 27 * i_irrp:27 * (i_irrp + 1)] += block

    return np.delete(B - B[i_conf_ref], i_conf_ref, axis=0).reshape(-1, 27 * len(irr_pairs))


def solve(A: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Least squares; minimum-norm solution with a warning if the data do not determine every element."""
    x, _, rank, _ = np.linalg.lstsq(A, y, rcond=None)
    if rank < A.shape[1]:
        print(f"Warning: {A.shape[1] - rank} combination(s) of the {A.shape[1]} unknowns "
              f"are not determined by the data; returning the minimum-norm solution.")
    return x


def check_fitting_quality(x: np.ndarray, y: np.ndarray, A: np.ndarray):
    """Write deviation statistics to fit.log and plot DFT vs fit."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y_fit = A @ x
    dy = np.abs(y_fit - y)

    # Percentages: the first element is always included, the rest only where y >= 0.001
    with np.errstate(divide="ignore", invalid="ignore"):
        dy_per = np.concatenate([[100 * dy[0] / y[0]], 100 * dy[1:][y[1:] >= 0.001] / y[1:][y[1:] >= 0.001]])

    with open("fit.log", "w") as f:
        f.write(f"Maximal deviation: {dy.max():12.6f} e Ang.\n")
        f.write(f"Average deviation: {dy.mean():12.6f} e Ang.\n")
        f.write(f" std of deviation: {dy.std(ddof=1):12.6f} e Ang.\n")
        f.write("\n")
        f.write(f"Maximal deviation in percentage: {dy_per.max():12.6f} \n")
        f.write(f"Average deviation in percentage: {dy_per.mean():12.6f} \n")
        f.write(f" std of deviation in percentage: {dy_per.std(ddof=1):12.6f} \n")

    index = np.arange(1, len(y) + 1)
    fig, ax = plt.subplots()
    ax.scatter(index, y, label="DFT")
    ax.scatter(index, y_fit, label="Fit")
    ax.set_xlabel(r"Index of $\delta P_\alpha$ for all spin configurations")
    ax.set_ylabel(r"$\delta P_\alpha$ (e $\AA$)")
    ax.legend(loc="lower right", frameon=False)
    fig.savefig("fitting_quality.pdf")
    plt.close(fig)


def save_P(x: np.ndarray, inp: dict, fname: str):
    """Write the 3x3 matrix P_gamma for every gamma and irreducible pair."""
    with open(fname, "w") as f:
        for i_irrp, k_pair in enumerate(inp["irr_pairs"]):
            i_site, j_site = inp["pairs"][k_pair - 1]
            f.write(f"Pair {i_site:3d} {j_site:3d}\n\n")
            for gamma, direction in enumerate("xyz"):
                f.write(f"The {direction} component of coupling matrix\n\n")
                for alpha in range(3):
                    index = 27 * i_irrp + 9 * gamma + 3 * alpha
                    f.write("{:12.6f} {:12.6f} {:12.6f}\n".format(*x[index:index + 3]))
                f.write("\n")
            f.write("\n")


def decompose_P(x: np.ndarray):
    """Split P into isotropic (J), antisymmetric (D), and anisotropic symmetric (Gamma) parts.

    P^J_{gamma alpha alpha} = trace(P_gamma) / 3, P^D = (P - P^T) / 2 over (alpha, beta),
    P^Gamma = P - P^J - P^D.
    """
    P = x.reshape(-1, 3, 3, 3)
    x_J = np.einsum("kg,ab->kgab", np.trace(P, axis1=2, axis2=3) / 3, np.eye(3)).flatten()
    x_D = ((P - P.transpose(0, 1, 3, 2)) / 2).flatten()
    x_Gamma = x - x_J - x_D
    return x_J, x_D, x_Gamma


def get_pi_and_d(x_J: np.ndarray, x_D: np.ndarray, inp: dict, fname: str):
    """Write pi (S_i . S_j term) and d_gamma (S_i x S_j term) for every irreducible pair.

    pi_gamma = P^J_{gamma xx}; d_{gamma x} = P^D_{gamma yz}, d_{gamma y} = P^D_{gamma zx}, d_{gamma z} = P^D_{gamma xy}.
    """
    P_J = x_J.reshape(-1, 3, 3, 3)
    P_D = x_D.reshape(-1, 3, 3, 3)
    with open(fname, "w") as f:
        for i_irrp, k_pair in enumerate(inp["irr_pairs"]):
            i_site, j_site = inp["pairs"][k_pair - 1]
            f.write(f"Pair {i_site:3d} {j_site:3d}\n\n")
            f.write("pi = {:12.6f} {:12.6f} {:12.6f}\n\n".format(*P_J[i_irrp, :, 0, 0]))
            for gamma, direction in enumerate("xyz"):
                d = P_D[i_irrp, gamma]
                f.write(f"d{direction} = {d[1, 2]:12.6f} {d[2, 0]:12.6f} {d[0, 1]:12.6f}\n")
            f.write("\n\n")


def get_P_in_local_reference_frames(x: np.ndarray, inp: dict) -> np.ndarray:
    """x_loc_{gamma alpha beta} = sum rmat_{gamma gamma'} rmat_{alpha alpha'} rmat_{beta beta'} x_{gamma' alpha' beta'}."""
    P = x.reshape(-1, 3, 3, 3)
    P_loc = np.array([np.einsum("gh,ai,bj,hij->gab", r, r, r, P[k]) for k, r in enumerate(inp["rmats"])])
    return P_loc.flatten()


def save_results(x: np.ndarray, inp: dict, modifier: str):
    save_P(x, inp, "P_total" + modifier + ".log")
    x_J, x_D, x_Gamma = decompose_P(x)
    save_P(x_J, inp, "P_J" + modifier + ".log")
    save_P(x_D, inp, "P_D" + modifier + ".log")
    save_P(x_Gamma, inp, "P_Gamma" + modifier + ".log")
    get_pi_and_d(x_J, x_D, inp, "pi_and_d" + modifier + ".log")


def print_relative_contribution(y_per: np.ndarray, f, tag: str):
    f.write(f"Maximal contribution from {tag:>8s}: {np.nanmax(y_per):6.1f} percent.\n")
    f.write(f"Average contribution from {tag:>8s}: {np.nanmean(y_per):6.1f} percent.\n")
    f.write(f" std of contribution from {tag:>8s}: {np.nanstd(y_per, ddof=1):6.1f} percent.\n")
    f.write("\n")


def get_relative_contribution(x: np.ndarray, A: np.ndarray, inp: dict):
    """Write the relative contributions of each irreducible pair and of its J, D, Gamma parts."""
    n_irrp = len(inp["irr_pairs"])
    blocks = [slice(27 * k, 27 * (k + 1)) for k in range(n_irrp)]
    y_irr = [np.abs(A[:, b] @ x[b]) for b in blocks]
    y_sum_irr = np.sum(y_irr, axis=0)
    x_J, x_D, x_Gamma = decompose_P(x)

    with np.errstate(divide="ignore", invalid="ignore"), open("relative_contribution.log", "w") as f:
        for k, b in enumerate(blocks):
            i_site, j_site = inp["pairs"][inp["irr_pairs"][k] - 1]
            f.write(f"Pair {i_site:3d} {j_site:3d}: relative contribution = "
                    f"{np.nanmean(100 * y_irr[k] / y_sum_irr):6.1f} percent.\n\n")
            y_J = np.abs(A[:, b] @ x_J[b])
            y_D = np.abs(A[:, b] @ x_D[b])
            y_Gamma = np.abs(A[:, b] @ x_Gamma[b])
            y_sum = y_J + y_D + y_Gamma
            print_relative_contribution(100 * y_J / y_sum, f, "P_J")
            print_relative_contribution(100 * y_D / y_sum, f, "P_D")
            print_relative_contribution(100 * y_Gamma / y_sum, f, "P_Gamma")


def save_P_for_one_pair_for_latex(f, x: np.ndarray, i_irrp: int):
    x = 1e5 * x
    for gamma, direction in enumerate("xyz"):
        f.write(f"The {direction} component of coupling matrix\n\n")
        for alpha in range(3):
            index = 27 * i_irrp + 9 * gamma + 3 * alpha
            end = " \\\\\n" if alpha != 2 else "\n"
            f.write("{:6.1f} & {:6.1f} &  {:6.1f}".format(*x[index:index + 3]) + end)
        f.write("\n")


def save_P_for_latex(x: np.ndarray, inp: dict):
    """Write P and its J, D, Gamma parts (in units of 1e-5 e Ang) as LaTeX table rows."""
    x_J, x_D, x_Gamma = decompose_P(x)
    with open("P_for_latex.log", "w") as f:
        for i_irrp, k_pair in enumerate(inp["irr_pairs"]):
            i_site, j_site = inp["pairs"][k_pair - 1]
            f.write(f"Pair {i_site:3d} {j_site:3d}\n\n")
            for title, part in [("Total ", x), ("The isotropic term", x_J),
                                ("The antisymmetric term", x_D), ("The anisotropic symmetric term", x_Gamma)]:
                f.write(f"{title}\n\n")
                save_P_for_one_pair_for_latex(f, part, i_irrp)
                f.write("\n")


def main():
    inp = read_yaml()
    S, P = read_dipoles(inp["spins"])

    y = get_y(P)
    A = get_A(S, inp)
    x = solve(A, y)

    # check_fitting_quality(x, y, A)
    # save_results(x, inp, "")

    x_loc = get_P_in_local_reference_frames(x, inp)

    # save_results(x_loc, inp, "_loc")
    # get_relative_contribution(x, A, inp)

    save_P_for_latex(x_loc, inp)
    return x, x_loc


if __name__ == "__main__":
    if len(sys.argv) > 1:
        os.chdir(sys.argv[1])
    main()
