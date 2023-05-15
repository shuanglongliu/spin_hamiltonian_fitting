import numpy as np

if __name__ == "__main__":

    data = np.loadtxt("energies.dat")
    m, n = data.shape
    l = m // 2
    data[0:l, n-1] = (data[0:l, n-1] + data[l:m, n-1])/2
    np.savetxt("energies_avergies.dat", data[0:l, :], fmt="%15.8f")

