"""
Calculate exchange coupling constants from DFT energies of different spin configurations.

Spin Hamiltonian: H = factor * sum_k Jk * sum_{ij in k} Si.Sj

For each magnetic configuration and each Jk:
    sumSS_k = sum_{ij in k} Si.Sj

Solves the linear system: factor * C * J = dE
where C is the matrix of sumSS differences and dE is the vector of energy differences.
"""

import numpy as np
import yaml


def read_data() -> np.ndarray:
    """Read magnetic configurations and energies from energies.dat.

    The first line must contain the energy unit ("eV" or "Ha").
    Energies are offset by the first entry and converted to cm^-1.

    Returns
    -------
    np.ndarray
        Array where each row is a magnetic configuration. The spin values occupy
        the leading columns and the last column holds the relative energy in cm^-1.
    """
    with open("energies.dat") as f:
        lines = f.readlines()

    if "eV" in lines[0]:
        unit_factor = 1000 * 8.06554393738
    elif "Ha" in lines[0]:
        unit_factor = 27.211396641308 * 1000 * 8.06554393738
    else:
        raise ValueError("Unknown energy unit in energies.dat (expected 'eV' or 'Ha').")

    data = np.loadtxt("energies.dat")
    data[:, -1] = (data[:, -1] - data[0, -1]) * unit_factor
    return data


def read_pairs() -> dict:
    """Read spin pairs for each J from pairs.yaml.

    Returns
    -------
    dict
        Mapping of J labels (e.g. "J1", "J2") to lists of 1-indexed spin-index pairs.
    """
    with open("pairs.yaml") as f:
        return yaml.safe_load(f)


def get_sumSS_vector(data: np.ndarray, pairs: list) -> np.ndarray:
    """Compute sum_{ij} Si.Sj for one J over all magnetic configurations.

    Parameters
    ----------
    data : np.ndarray
        Configuration array; rows are configurations, spin values in all but last column.
    pairs : list of [int, int]
        1-indexed spin-index pairs that share this J.

    Returns
    -------
    np.ndarray
        sumSS[k] = sum of Si.Sj products for configuration k.
    """
    sumSS = np.zeros(data.shape[0])
    for i, j in pairs:
        sumSS += data[:, i - 1] * data[:, j - 1]
    return sumSS


def get_sumSS_matrix(data: np.ndarray, pairs_dict: dict) -> np.ndarray:
    """Build the sumSS matrix for all Js and all magnetic configurations.

    Parameters
    ----------
    data : np.ndarray
        Configuration array.
    pairs_dict : dict
        Mapping of J labels to spin-index pair lists.

    Returns
    -------
    np.ndarray
        Shape (n_J, n_config). Entry [i, j] is sum_{kl} Sk.Sl for Ji
        evaluated on configuration j.
    """
    J_labels = list(pairs_dict.keys())
    sumSS_matrix = np.zeros((len(J_labels), data.shape[0]))
    for i, label in enumerate(J_labels):
        sumSS_matrix[i] = get_sumSS_vector(data, pairs_dict[label])
    return sumSS_matrix


def get_Js(data: np.ndarray, sumSS_matrix: np.ndarray, factor: float) -> np.ndarray:
    """Solve for exchange coupling constants from the linear system factor * C * J = dE.

    Energy differences are taken with respect to the last magnetic configuration.
    An exact solve is used when the system is square; least-squares when overdetermined.

    Parameters
    ----------
    data : np.ndarray
        Configuration array; last column holds relative energies in cm^-1.
    sumSS_matrix : np.ndarray
        Shape (n_J, n_config) matrix from get_sumSS_matrix.
    factor : float
        Prefactor in the spin Hamiltonian (e.g. -2 for H = -2J Si.Sj).

    Returns
    -------
    np.ndarray
        Exchange coupling constants in cm^-1.

    Raises
    ------
    ValueError
        If the system is underdetermined (fewer equations than unknowns).
    """
    n_ref = data.shape[0] - 1   # index of reference configuration
    n_eq = n_ref                 # number of equations after differencing
    n_J = sumSS_matrix.shape[0]

    if n_eq < n_J:
        raise ValueError(
            f"Underdetermined system: {n_eq} equation(s) for {n_J} unknown(s)."
        )

    # C[i, j] = sumSS_matrix[j, i] - sumSS_matrix[j, n_ref]
    C = (sumSS_matrix[:, :n_eq] - sumSS_matrix[:, [n_ref]]).T  # shape (n_eq, n_J)
    dE = data[:n_eq, -1] - data[n_ref, -1]

    if n_eq == n_J:
        Js = np.linalg.solve(C, dE) / factor
    else:
        Js = np.linalg.lstsq(C, dE, rcond=None)[0] / factor

    values = "  ".join(f"{J:6.2f}" for J in Js)
    print(f"Js = [{values}] cm^-1")
    return Js


if __name__ == "__main__":
    data = read_data()
    pairs_dict = read_pairs()
    sumSS_matrix = get_sumSS_matrix(data, pairs_dict)
    Js = get_Js(data, sumSS_matrix, factor=-2)
