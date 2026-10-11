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


def read_data(energy_threshold=None) -> np.ndarray:
    """Read magnetic configurations and energies from energies.dat.

    The first line must contain the energy unit ("eV" or "Ha").
    Energies are offset by the first entry and converted to cm^-1.

    Input
    -----
    energy_threshold : float
        Optional threshold in cm^-1 to exclude high-energy outliers relative to the lowest-energy configuration

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

    # Auto detect and remove high-energy outliers relative to the lowest-energy configuration 
    if energy_threshold is not None:
        min_energy = np.min(data[:, -1])
        energy_diff = data[:, -1] - min_energy
        mask = energy_diff <= energy_threshold
        if not np.all(mask):
            print(f"Excluding {np.sum(~mask)} configuration(s) with energy > {energy_threshold} cm^-1 above the minimum.")
            data = data[mask]
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
        # Report m, and n for exact solve
        print(f"Square system: {n_eq} equations for {n_J} unknowns. Performing exact solve.")
        Js = np.linalg.solve(C, dE) / factor
    else:
        # Report m, and n for least-squares fit
        print(f"Overdetermined system: {n_eq} equations for {n_J} unknowns. Performing least-squares fit.")
        Js = np.linalg.lstsq(C, dE, rcond=None)[0] / factor
        # Report fitting quality
        dE_fit = factor * C @ Js
        residual = np.linalg.norm(dE - dE_fit)
        print(f"Least-squares fit: residual = {residual:.2e} cm^-1")
        # Plot fitted vs actual energies
        import matplotlib.pyplot as plt
        plt.scatter(dE, dE_fit, label="Data vs Fit")
        plt.plot(dE, dE, "r--", label="Ideal Fit")
        plt.xlabel("Actual dE (cm^-1)")
        plt.ylabel("Fitted dE (cm^-1)")
        plt.legend()
        plt.title("Actual vs Fitted Energy Differences")
        plt.savefig("fit_quality.png")
        print("Fit quality plot saved to fit_quality.png.")

    # Save Js to a csv file
    with open("Js.csv", "w") as f:
        f.write(f"# Js in cm^-1, factor = {factor}\n")
        f.write(",".join(f"J{i+1}" for i in range(len(Js))) + "\n")
        f.write(",".join(f"{J:.2f}" for J in Js) + "\n")
    print(f"Calculated Js saved to Js.csv with factor = {factor}.")


if __name__ == "__main__":
    data = read_data(energy_threshold=1000)
    pairs_dict = read_pairs()
    sumSS_matrix = get_sumSS_matrix(data, pairs_dict)
    Js = get_Js(data, sumSS_matrix, factor=-2)
