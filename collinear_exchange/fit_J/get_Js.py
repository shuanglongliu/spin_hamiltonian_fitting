# Calculate the exchange coupling constants from DFT energies of different spin configurations.
# Spin Hamiltonian: H = factor*J1 Si.Sj + factor*J2 Si.Sj + ...
# sumSS1 = sum Si.Sj for J1
# sumSS2 = sum Si.Sj for Jk
# ... for all Js. 
# Then solve the equation factor Cs Js = dEs, where Cs is the matrix of sumSSs and dEs is the vector of energy differences from E0.

import numpy as np

def read_data():
    data = np.loadtxt("energies.dat")
    return data

def get_sumSS(data, pairs):
    # H = factor*J1 sum_ij Si.Sj + factor*J2 sum_ij Si.Sj + ...
    # pairs = [[0,1], [0,2], ...] for a certain J.
    # sumSS = sum_ij Si.Sj for the corresponding J. 
    # One sumSS for each magnetic configuration (each row of data).
    npair = len(pairs)
    As = []
    for i in range(data.shape[0]):
        A = 0
        for ipair in range(npair):
            A = A + data[i, pairs[ipair][0]] * data[i, pairs[ipair][1]]
        As.append(A)
    return As

def get_Js(data, sumSSs, factor):
    # H = factor*J1 sum_ij Si.Sj + factor*J2 sum_ij Si.Sj + ...
    # Set up and solve the equation factor Cs Js = dEs
    # Cs[i, j] = sumSSs[j][n] - sumSSs[j][i]
    # sumSSs[i][m] the coefficient of Ji for the magnetic configuration m.
    n = data.shape[0] - 1 # number of linear equations after subtracting E0
    Cs = np.zeros((n, n))
    dEs = []
    for i in range(n):
        for j in range(n):
            Cs[i, j] = sumSSs[j][n] - sumSSs[j][i]
        dEs.append(data[n, -1] - data[i, -1])
    if n == len(sumSSs):
        Js =  8.06554393738 * 1000 * np.matmul(np.linalg.inv(Cs), dEs)/factor
    elif n > len(sumSSs):
        Js =  8.06554393738 * 1000 * np.linalg.lstsq(Cs, dEs, rcond=None)[0]/factor
    else:
        raise ValueError("Number of equations must be greater than or equal to the number of Js.")
    print(("Js = [" +  n*" {:6.2f} " + "] cm^-1").format(*Js, *Js))
    return Js

if __name__ == "__main__":

    data = read_data()
    sumSS1 = get_sumSS(data, [[0,1]])
    sumSS2 = get_sumSS(data, [[0,2]])
    sumSS3 = get_sumSS(data, [[0,3]])
    sumSS4 = get_sumSS(data, [[1,2]])
    sumSS5 = get_sumSS(data, [[1,3]])
    sumSS6 = get_sumSS(data, [[2,3]])
    sumSS7 = get_sumSS(data, [[0,4]])
    sumSS8 = get_sumSS(data, [[0,5]])
    sumSS9 = get_sumSS(data, [[0,6]])
    sumSS10 = get_sumSS(data, [[4,5]])
    sumSS11 = get_sumSS(data, [[4,6]])
    sumSS12 = get_sumSS(data, [[5,6]])
    Js = get_Js(data, [sumSS1, sumSS2, sumSS3, sumSS4, sumSS5, sumSS6, sumSS7, sumSS8, sumSS9, sumSS10, sumSS11, sumSS12], -2)

